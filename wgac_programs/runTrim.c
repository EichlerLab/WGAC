#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#define FALSE 0
#define TRUE !FALSE
#define WORKTAG 1
#define DIETAG 2
#define _NULL 0
#define BIG_STRSIZE 4096
#define MED_STRSIZE 1024
#define SMALL_STRSIZE 256


struct _fileliststruct{

	char *name;
	int started;
      int finished;
	int count;
};

typedef struct _fileliststruct filelist;

/* MPI FUNCTIONS */
static void master(char *,char *,char *,char *,char *);
static void slave();
char *get_next_work_item(char *,char *,char *,char *,char *,filelist **);
/* MPI FUNCTIONS */


//Functions related to file processing
int file_select();
extern  int alphasort();
filelist *processtheinputlist_and_initialize(char *inputdir);
//Functions related to file processing



int main(int argc, char **argv)
{

  char *parameters;
  char *inputdir;
  char *fastadir;
  char *outputdir;
  char *jobid;
  int myrank;
  
  parameters = (char *) malloc(sizeof(char)*BIG_STRSIZE);
  inputdir = (char *) malloc(sizeof(char)*BIG_STRSIZE);
  fastadir = (char *) malloc(sizeof(char)*BIG_STRSIZE);
  outputdir = (char *) malloc(sizeof(char)*BIG_STRSIZE);
  jobid = (char *) malloc(sizeof(char)*SMALL_STRSIZE);
 
  strcpy(parameters,argv[1]);
  strcpy(inputdir,argv[2]);
  strcpy(fastadir,argv[3]);
  strcpy(outputdir,argv[4]);
  strcpy(jobid,argv[5]);
  
  printf("SABA: parameters:%s\n inputdir:%s\n fastadir:%s\n outputdir:%s\n",parameters,inputdir,fastadir,outputdir);

 /* Initialize MPI */

  MPI_Init(&argc,&argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    master(inputdir,fastadir,outputdir,parameters,jobid);
  } else {

    printf("slave is started\n");
    slave();
  }
                                                               
  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


static void master(char *inputdir,char *fastadir,char *outputdir,char *parameters, char *jobid)
{
  int ntasks, rank;
  char *result;
  MPI_Status status;
  
  filelist *myfilelist;
  char *work;
  char *diestring;
  DIR *dfd;

  int count;
  struct direct **files;
  /* Find out how many processes there are in the default
     communicator */
 diestring=(char *)malloc(sizeof(char)*SMALL_STRSIZE);
 result=(char *)malloc(sizeof(char)*MED_STRSIZE);
 work=(char *)malloc(sizeof(char)*BIG_STRSIZE);
 MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  /* Seed the slaves; send one unit of work to each slave. */
	strcpy(diestring,"XXX"); //if work is null, it means the number of jobs are more than the processors 


  printf("MASTER Master received the inputs\ninputdir is %s \nparameters are %s\n",inputdir,parameters);
  

  if ((dfd=opendir(inputdir))==NULL)
  {
        printf("Can't access %s \nSo I am quitting\n",inputdir);
        exit(0); //

  }
  closedir(dfd);

  count=scandir(inputdir,&files,file_select,alphasort);

  myfilelist=(filelist *) malloc(sizeof(filelist)*(count+1));
  myfilelist=processtheinputlist_and_initialize(inputdir);


  for (rank = 1; rank < ntasks; ++rank) {

    /* Find the next item of work to do */

    printf("MASTER Forloop ranks\n");  
    work = get_next_work_item(inputdir,fastadir,outputdir,parameters,jobid,&myfilelist);
    /* Send it to each rank */
    if (work == _NULL)
    { 
        MPI_Send(diestring,             /* message buffer */
                BIG_STRSIZE,                 /* one data item */
                MPI_CHAR,           /* data item is an integer */
                rank,              /* destination process rank */
                DIETAG,           /* user chosen message tag */
                MPI_COMM_WORLD);   /* default communicator */
        printf("DIEstring sent \n");



    }
    else 
    {
    	printf("MASTER work issssssssssssss %s \n",work);
    	MPI_Send(work,             /* message buffer */
       		BIG_STRSIZE,                 /* one data item */
             	MPI_CHAR,           /* data item is an integer */
             	rank,              /* destination process rank */
             	WORKTAG,           /* user chosen message tag */
             	MPI_COMM_WORLD);   /* default communicator */
     	strcpy(work,"\0"); //initialize it
     	printf("MASTER sent success \n");
    }
    fflush(stdout);

  }

  printf("MASTER First phase is completed\n");
  /* Loop over getting new work requests until there is no more work
     to be done */
  
  work = get_next_work_item(inputdir,fastadir,outputdir,parameters,jobid,&myfilelist);

  while (work != _NULL) {

    /* Receive results from a slave */

    MPI_Recv(result,           /* message buffer */
             MED_STRSIZE,                 /* one data item */
             MPI_CHAR,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

    //process the result, update &myfilelist
    /* Send the slave a new work unit */

    MPI_Send(work,             /* message buffer */
             BIG_STRSIZE,                 /* one data item */
             MPI_CHAR,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    /* Get the next unit of work to be done */
    work = get_next_work_item(inputdir,fastadir,outputdir,parameters,jobid,&myfilelist);
    fflush(stdout);

    
  }
  printf("MASTER While is skipped\n");
  printf("still here\n");
  /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */

  for (rank = 1; rank < ntasks; ++rank) {
    printf("Rank %d is done\n",rank);	
    MPI_Recv(result, MED_STRSIZE, MPI_CHAR, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  	 fflush(stdout);

  }
  printf("MASTER All the results are collected, now sending die signals to processors\n");
  /* Tell all the slaves to exit by sending an empty message with the
     DIETAG. */
  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Send(diestring, BIG_STRSIZE, MPI_CHAR, rank, DIETAG, MPI_COMM_WORLD);
  }
 printf("MASTER Program finished succesfully \n");
}


static void slave()
{
  char  work[BIG_STRSIZE];
  MPI_Status status;
  char *command=(char *)malloc(sizeof(char)*BIG_STRSIZE);
  char *cdcommand=(char *)malloc(sizeof(char)*BIG_STRSIZE);
  char *result=(char *)malloc(sizeof(char)*MED_STRSIZE);
  //char *outdir=(char *)malloc(sizeof(char)*MED_STRSIZE);
  char *mvop=(char *)malloc(sizeof(char)*MED_STRSIZE);
  char *myjobid=(char *)malloc(sizeof(char)*MED_STRSIZE);
  char *currentfile=(char *)malloc(sizeof(char)*MED_STRSIZE);
  char fromfile[1000];
  char tofile[1000];

  sprintf(result,"%s","noworkassigned");
  while (1) {

    /* Receive a message from the master */
    MPI_Recv(work, BIG_STRSIZE, MPI_CHAR, 0, MPI_ANY_TAG,
	             MPI_COMM_WORLD, &status);

    /* Check the tag of the received message. */

   printf("\n\n\nSLAVE STARTED \n work is  %s\n",work);
  //either master send me a "DIE" command, or it is the beginning and number of works are less than number of processors
  //no work to be done so I just return 
   if ((status.MPI_TAG == DIETAG) ) {
     printf("DIE SIGNAL SENT TO ME\n"); 
     MPI_Send(result, MED_STRSIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

     return;
    }
   if (!strcmp(work,"XXX"))
    {
	printf("No more work to be done so quitting...\n");
	MPI_Send(result, MED_STRSIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	return;

    }


    /* Do the work */
    command=strtok(work,"\n");
    myjobid=strtok(_NULL,"\n");
    currentfile=strtok(_NULL,"\n");

    system(command); //DO THE MRSFAST HERE
    
    printf("Work is %s \n",work);
    printf("SLAVE FINISHED, Sending the result \n");
    /* Send the result back */
    MPI_Send(result, MED_STRSIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    printf("Successfully sent \n");

    //free(command);
    //free(cdcommand);
    //free(result);
    //free(outdir);
  }
}


char *get_next_work_item(char *inputdir,char *fastadir,char *outputdir,char *parameters,char *jobid,filelist **myfilelist)
{
  /* Fill in with whatever is relevant to obtain a new unit of work
     suitable to be given to a slave. */
//Input: filedir
//Input: outdir

	char *currentfile;
  	char *out;
  	int i,count;
  	char *myjobid;

	currentfile=(char *)malloc(sizeof(char)*MED_STRSIZE);
	out=(char *)malloc(sizeof(char)*BIG_STRSIZE);
	myjobid=(char *)malloc(sizeof(char)*SMALL_STRSIZE);
	sprintf(myjobid,"%s.1.all.q",jobid);

	count = (*myfilelist)[0].count;

 	for (i=0;i<count;i++) 
 	{

		if (((*myfilelist)[i].finished==0) && ((*myfilelist)[i].started==0)) //if that file hasn't been touched by other processor yet
	 	{

			strcpy(currentfile,(*myfilelist)[i].name); //get the file name
			sprintf(out,"%s %s/%s %s %s\n%s\n%s\n",parameters,inputdir,currentfile,fastadir,outputdir,myjobid,currentfile);
			printf("SABA: out:%s",out);

			(*myfilelist)[i].started=1; //mark that file as started.
			return out;
        	}
	}

	//free(currentfile);
	//free(out);
	//free(myjobid);
	return _NULL; //no files to go but how about the ones that are started but not finished
}
filelist *processtheinputlist_and_initialize(char *inputdir)
{ 
//scandir routine is from www.cs.cf.ac.uk/Dave/C/node20.html
  int count,i;
  struct direct **files;
  filelist *myfilelist;
  char *outfile;
  DIR *dfd;
  FILE *outputfilepointer;
  outfile=(char *)malloc(sizeof(char)*MED_STRSIZE);   
    
    
    //check if the dir exists??? ADD MORE ERROR CHECKING HERE!!!!!!
  printf("Current Working Directory = %s\n",inputdir);
        
         
  if ((dfd=opendir(inputdir))==NULL)
  {
        printf("Can't access %s \nSo I am quitting\n",inputdir);
        exit(0); //
                
  }
  closedir(dfd);

  count = scandir(inputdir, &files, file_select, alphasort);


  printf("Number of files in the list are  %d \n",count);
  myfilelist=(filelist *) malloc(sizeof(filelist)*(count+1));
  
 
 for (i=0;i<count;i++)
 {
     myfilelist[i].name=(char *)malloc(sizeof(char)*MED_STRSIZE);
     strcpy(myfilelist[i].name,files[i]->d_name);
     //sprintf(outfile,"%s/%s.out",outputdir,myfilelist[i].name);  

     //printf("outfile is %s\n",outfile);

     if (strstr(files[i]->d_name, "summary.txt")){
       printf("PRE PROCESSING: this file  %s  a summary file.\n",files[i]->d_name);
       myfilelist[i].finished=1;
       myfilelist[i].started=1;
       myfilelist[i].count=count;
       continue;
     }


     if ((outputfilepointer=fopen(outfile,"r"))==NULL)
     {
        printf("PRE PROCESSING : this file  %s  hasn't processed yet\n",files[i]->d_name);
        myfilelist[i].finished=0;
        myfilelist[i].started=0;
        myfilelist[i].count=count;
       
     }
     else
     { 
	fclose(outputfilepointer);
        printf("PRE PROCESSING: this file  %s  already processed\n",files[i]->d_name);
        myfilelist[i].finished=1;
        myfilelist[i].started=1;
	myfilelist[i].count=count;

     }  
                
 }  


 fflush(stdout);

 //free(outfile);  
 return myfilelist;
  
} 
int file_select(struct direct *entry)
                                                                                                                                                             
{
        if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0))
                return (FALSE);
        else
                return (TRUE);
}






