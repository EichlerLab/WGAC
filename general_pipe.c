#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>


#define FALSE 0
#define TRUE !FALSE
#define WORKTAG 1
#define DIETAG 2
#define _NULL 0


enum numbers { HUGESTR = 2048, BIGSTR = 1024, MEDSTR = 512, SMLSTR = 128 };

//To compile:mpicc -g -o general_pipe general_pipe.c


struct _fileliststruct
{
  char *name; //name of the file
  int started; //0 means file haven't touched yet,1 means either it is still processing or finished
  int finished; //1 means file has finished, 0 means it is either still processing or not touched yet
};


typedef struct _fileliststruct filelist;


///GLOBAL VARIABLES
//These are stable anyway, so there is no harm in defining them global
int  INPUTSIZECOUNT;
char **handy_param;     //stores char *s of parameters seperated by dummy_in and dummy_out
int  cmd_len = 0;       //# of unit separated by space in command line
int  dbug = 1;
char *SEPERATOR = "******";  //used to seperate output filename from the command to run in slave


/* MPI FUNCTIONS */
static void master(char *, char *, char *, char *, char *);
static void slave();
char *get_next_work_item(char *, filelist **, char *, char *, char *);
/* MPI FUNCTIONS   */



//Functions related to file processing
int file_select();
extern  int alphasort();
filelist *processtheinputlist_and_initialize(char *inputdir, char *outdir, char *suffix);
void process_cmd(char *param, char *indir, char *outdir);


int main(int argc, char **argv)
{ 
  if( argc < 6 )
    {
      printf("%s usage     parameters     inputdir     outputdir   suffix   tmpdir\nYou only give %d args\n", argv[0], argc);
      int i = 0;
      for(; i < argc; i++)
	printf("arg %d: '%s'\n", i, argv[i]);
      exit(0);
    }
  char *parameters;
  char *inputdir;
  char *outputdir;
  char *suffix;           //suffix following input file name to give output filename, if using dot, include dot
  char *tmpdir;           //tmp local dir to put files before write to real out dir to speed up gpfs ops.
  char **MPIargv;
  int  myrank;
  parameters = (char *) malloc(sizeof(char)* BIGSTR);
  inputdir   = (char *) malloc(sizeof(char)* BIGSTR);
  outputdir  = (char *) malloc(sizeof(char)* BIGSTR);
  suffix     = (char *) malloc(sizeof(char)* SMLSTR);
  tmpdir     = (char *) malloc(sizeof(char)* BIGSTR);
  // processing command, remember strtok chop tokens off the src, so can't be processed >1
  strcpy(parameters, argv[1]);
  strcpy(inputdir,   argv[2]);
  strcpy(outputdir,  argv[3]);
  strcpy(suffix,     argv[4]);
  strcpy(tmpdir,     argv[5]);
  /* Initialize MPI */
  MPI_Init(&argc,&argv);
  /* Find out my identity in the default communicator */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
    master(inputdir, outputdir, parameters, suffix, tmpdir);
  else
    {
      printf("slave is started\n");
      slave();
    }
  /* Shut down MPI */
  MPI_Finalize();
  return 0;
}


static void master(char *inputdir, char *outputdir, char *parameters, char *suffix, char *tmpdir)
{
  int ntasks, rank;
  char *result;
  MPI_Status status;
  filelist *myfilelist;
  char *work;
  char *diestring;
  char *command;
  DIR *dirpointer1,*dirpointer2;
  /* Find out how many processes there are in the default communicator */
  diestring = (char *)malloc(sizeof(char)*MEDSTR);
  result    = (char *)malloc(sizeof(char)*MEDSTR);
  work      = (char *)malloc(sizeof(char)*HUGESTR);
  command   = (char *)malloc(sizeof(char)*HUGESTR);
  //process command line to initialize cmd unit array and remove trailing backslash of dirs
  process_cmd(parameters, inputdir, outputdir);
  //try input dir and output dir
  if (((dirpointer1 = opendir(inputdir)) == NULL) || ((dirpointer2 = opendir(outputdir)) == NULL) )
    {
      printf("At least one of the input directories given as input was not accessible\n");
      exit(0);
    }
  closedir(dirpointer1);
  closedir(dirpointer2);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  strcpy(diestring,"XXX");
  //if work is null, it means the number of jobs are more than the processors
  //if there is a output file in the output dir, not included in the to be worked on
  myfilelist = processtheinputlist_and_initialize(inputdir, outputdir, suffix);
  /* Seed the slaves; send one unit of work to each slave. */
  for (rank = 1; rank < ntasks; ++rank)
    {
      /* Find the next item of work to do */
      printf("MASTER FORLOOP RANK Number %d \n", rank);
      work = get_next_work_item(inputdir, &myfilelist, outputdir, suffix, tmpdir); 
      if (work == _NULL)
	{
	  MPI_Send(diestring,             /* message buffer */
		   HUGESTR,                 /* one data item */
		   MPI_CHAR,           /* data item is an integer */
		   rank,              /* destination process rank */
		   DIETAG,           /* user chosen message tag */
		   MPI_COMM_WORLD);   /* default communicator */
	  printf("No work to do, DIEstring sent \n");
	}
      else
	{
	  printf("MASTER work is %s \n",work);
	  MPI_Send(work,             /* message buffer */
		   HUGESTR,                 /* one data item */
		   MPI_CHAR,           /* data item is an integer */
		   rank,              /* destination process rank */
		   WORKTAG,           /* user chosen message tag */
		   MPI_COMM_WORLD);   /* default communicator */
	  strcpy(work, "\0"); //initialize it
	}
    }
  
  printf("MASTER First phase is completed\n");
  /* Loop over getting new work requests until there is no more work  to be done */
  

  work = get_next_work_item(inputdir, &myfilelist, outputdir, suffix, tmpdir);
  while (work != _NULL)
    {
      
      /* Receive results from a slave */
      MPI_Recv(result,           /* message buffer */
	       MEDSTR,                 /* one data item */
	       MPI_CHAR,        /* of type double real */
	       MPI_ANY_SOURCE,    /* receive from any sender */
	       MPI_ANY_TAG,       /* any type of message */
	       MPI_COMM_WORLD,    /* default communicator */
	       &status);          /* info about the received message */
      printf("%s\n", result);

      //process the result, update &myfilelist
      /* Send the slave a new work unit */
      
      MPI_Send(work,             /* message buffer */
	       HUGESTR,                 /* one data item */
	       MPI_CHAR,           /* data item is an integer */
	       status.MPI_SOURCE, /* to who we just received from */
	       WORKTAG,           /* user chosen message tag */
	       MPI_COMM_WORLD);   /* default communicator */
      
      /* Get the next unit of work to be done */
      work = get_next_work_item(inputdir, &myfilelist,  outputdir, suffix, tmpdir);
    }
  /* There's no more work to be done, so receive all the outstanding  results from the slaves. */



  for (rank = 1; rank < ntasks; ++rank)
    {
      printf("Rank %d is done\n",rank);
      MPI_Recv(result, MEDSTR, MPI_CHAR, MPI_ANY_SOURCE,  MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      printf("%s\n", result);
    }
  printf("MASTER All the results are collected, now sending die signals to processors\n");


  /* Tell all the slaves to exit by sending an empty message with the DIETAG. */
  for (rank = 1; rank < ntasks; ++rank)
    {
      MPI_Send(diestring, HUGESTR, MPI_CHAR, rank, DIETAG, MPI_COMM_WORLD);
    }
  printf("MASTER Program finished succesfully \n");
}



// let slave take arguments from main function does not work
static void slave()
{
  char work[HUGESTR];
  char *result;
  MPI_Status status;
  char *command, *newwork;
  int  rtn_state;
  char *currentfile, *outfile;
  FILE *temp;


  currentfile = (char *)malloc(sizeof(char)*HUGESTR);
  command     = (char *)malloc(sizeof(char)*HUGESTR);
  result      = (char *)malloc(sizeof(char)*MEDSTR);
  

  //owork=(char *)malloc(sizeof(char)*HUGESTR);
  while (1)
    {

      /* Receive a message from the master */
      MPI_Recv(work, HUGESTR, MPI_CHAR, 0, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &status);
      

      /* Check the tag of the received message. */
      printf("\n\n\nSLAVE STARTED \n work is  %s\n", work);
      //either master send me a "DIE" command, or it is the beginning and number of works are less than number of processors
      //no work to be done so I just return
      
      
      if ((status.MPI_TAG == DIETAG) )
	{
	  printf("DIE SIGNAL SENT TO ME\n");
	  return;
	}
      if (!strcmp(work, "XXX"))
	{
	  printf("No more work to be done so quitting...\n");
	  return;
	}
    

      //seperate the current output file name first by seperator
      if( (outfile = strstr(work, SEPERATOR)) == NULL )
	{
	  sprintf(result, "'%s' does not have '%s' as seperator to seperate output file and command. Not executed\n", work, SEPERATOR);
	  MPI_Send(result, MEDSTR, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
      bzero(command,   HUGESTR);
      strncpy(command, work, outfile - work);
      
      /* Do the work */
      newwork   = outfile + strlen(SEPERATOR);
      rtn_state = system(newwork);
      if( rtn_state != 0 )
	sprintf(result, "'%s' is not working right.      %s\n", newwork, strerror(errno));
      else
	{
	  rtn_state = system(command);
	  if( rtn_state != 0 )
	    sprintf(result, "mv outfile for work '%s' is not working right with the mv command as '%s'.     %s\n", newwork, command, strerror(errno));
	}
   

      //CHECK IF FILES EXIST HERE BEFORE MOVING@!!!!!!!!!!!!!
      printf(" SLAVE FINISHES \n\n");
      
      /* Send the result back */
      MPI_Send(result, MEDSTR, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
}






char *get_next_work_item(char *inputdir, filelist **myfilelist, char *outdir, char *suffix, char *tmpdir)
{
  char *currentfile;
  char *out;
  int i = 0;
  int j;
  int count;
  struct direct **files;
  DIR *dirpointer1;


  currentfile = (char *)malloc(sizeof(char)*MEDSTR);
  out         = (char *)malloc(sizeof(char)*HUGESTR);
  count       = INPUTSIZECOUNT;


  for (i=0; i<count; i++)
    {
      bzero(out, HUGESTR);

      
      if (  ( (*myfilelist)[i].finished == 0) && ( (*myfilelist)[i].started == 0)    ) //if that file hasn't been touched by other processor yet
	{
	  strcpy(currentfile,(*myfilelist)[i].name);                              //get the file name
	  //use '*' as delim to seperate the commands as * should not be in file name
	  sprintf(out,     "mv  %s/%s%s  %s%s", tmpdir, (*myfilelist)[i].name, suffix, outdir, SEPERATOR);	  
	  
	  //replace the dummy_in and dummy_out in the command with the real input and output
	  //this is only used for 1 input to 1 output
	  for(j = 0; j < cmd_len; j++)
	    {
	      if( strcmp(handy_param[j], "dummy_in" ) == 0 )
		{
		  strcat(out, inputdir);
		  strcat(out, "/");
		  strcat(out, currentfile);
		  strcat(out, " ");
		  continue;
		}
	      if( strcmp(handy_param[j], "dummy_out" ) == 0 )
		{
		  strcat(out, tmpdir);
		  strcat(out, "/");
		  strcat(out, currentfile);
		  strcat(out, suffix);
		  strcat(out, " ");
		  continue;
		}
	      
	      strcat(out, handy_param[j]);
	      strcat(out, " ");
	    }
	
	  //changed to add currentfile name at the end so 
	  (*myfilelist)[i].started = 1; //mark that file as started.
	  //!!!!!!: WHAT IF THAT FILE IS STARTED AND NEVER FINISHES STARTED FLAG IS GOING TO REMAIN 1,
	  //?? Right now only solution is to restart the rmasker_batch
	  printf("in get next item: %s\n", out); 
	  return out;
	}   
    }
  
  return _NULL; //no files to go but how about the ones that are started but not finished
}





filelist *processtheinputlist_and_initialize(char *inputdir, char *outdir, char *suffix)
{
  //scandir routine is from www.cs.cf.ac.uk/Dave/C/node20.html
  int count,i;
  struct direct **files;
  filelist *myfilelist;
  char     *outfile;
  DIR      *dfd;
  FILE     *tempfp;

  outfile = (char *) malloc(sizeof(char) * HUGESTR);


  //check if the dir exists??? ADD MORE ERROR CHECKING HERE!!!!!!
  if( dbug )
    printf("Current Working Directory = %s\n", inputdir);
  
  
  if ( (dfd=opendir(inputdir)) == NULL )
    {
      printf("Can't access %s \nSo I am quitting\n",inputdir);
      exit(0); //
    }
  count = scandir(inputdir, &files, file_select, alphasort);


  INPUTSIZECOUNT=count;
  printf("Number of files in the list are  %d \n",count);
  myfilelist=(filelist *) malloc(sizeof(filelist)*(count+1));


  for (i = 0; i < count; i++)
    {
      myfilelist[i].name = (char *)malloc(sizeof(char)*MEDSTR);
      strcpy(myfilelist[i].name, files[i]->d_name);
      bzero(outfile, HUGESTR);

      sprintf(outfile,  "%s/%s%s",      outdir,     myfilelist[i].name, suffix);    

      if ( (tempfp = fopen(outfile,"r")) == NULL )
	{
	  //if any of those files are missing, that means that repeatmasker part hasn't completed,so I initialize the flags;
	  printf("PRE PROCESSING : this file  %s  hasn't processed yet\n",files[i]->d_name);
   	  myfilelist[i].finished=0;
          myfilelist[i].started=0;
	}
      else //If all the output is there, it means I don't need to do some repeatmask anymore, so I flag them.
	{
	  fclose(tempfp);
	  printf("PRE PROCESSING: this file  %s  already processed\n", files[i]->d_name);
	  myfilelist[i].finished=1;
	  myfilelist[i].started=1;
	}
      
    }
  fflush(stdout);
  sleep(100);
  closedir(dfd);
  return myfilelist;
}


int file_select(struct direct *entry) 
{
  if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0))
    return (FALSE);
  else
    return (TRUE);
}



void process_cmd(char *param, char *indir, char *outdir)
{
  char p2[HUGESTR], parameters[HUGESTR];
  char *token;
  char see_in  = '0';    //if sees dummy_in  in cmd
  char see_out = '0';    //if sees dummy_out in cmd


  //get the number of tokens in cmd
  strcpy(p2, param);
  token = strtok(p2, " ");
  while( token != NULL )
    {
      cmd_len++;
      token = strtok(NULL, " ");
    }
  handy_param = (char **) malloc(sizeof(char *) * cmd_len);


  int i = 0;
  strcpy(parameters, param);
  token = strtok(parameters, " ");
  while( token != NULL )
    {
      char *cmd_unit   = (char *) malloc(sizeof(char)* BIGSTR);
      strcpy(cmd_unit, token);
      handy_param[i++] = cmd_unit;

      if( strcmp(cmd_unit, "dummy_in") == 0 )
	see_in  = '1';
      if( strcmp(cmd_unit, "dummy_out") == 0 )
	see_out = '1';
      
      token = strtok(NULL, " ");
    }


  if( dbug )
    {
      printf("output command unit  ---------------------------------------------\n");
      for(i = 0; i < cmd_len; i++)
	{
	  printf("%s\n", handy_param[i]);
	}
    }
  


  if( see_in == '0' || see_out == '0' )
    {
      printf("Please have dummy_in in the place of input file and dummy_out in the place of output file in the command line\n");
      exit(1);
    }

  
  //remove back slashes
  if( indir[ strlen(indir) -1 ] == '/' )
    indir[ strlen(indir) -1 ] = '\0';
  
  if( outdir[ strlen(outdir) -1 ] == '/' )
    outdir[ strlen(outdir) -1 ] = '\0';
  
}






/* MPI program notes:
   1. Master and Slave only share the CODE, not the process of the running program.
   So DO NOT expect anything global initialized by master can be accessable by slave
   The only way for them to communicate is thru the string passed

   2. Only use global variables when they : needed by all master used functions
                                            needed by all slave used functions
*/
