#include    <stdio.h>
#include    <ctype.h>
#include    <stdbool.h>
#include <locale.h>


int main( int argc, char* argv[] ) {
   char szLastHeader[10000];
   int nHeaderLength = 0;
   int nNumberOfBases = 0;
   bool bInHeader = false;
   long long llNumberOfBasesAllSequences = 0;
   char c = 'x';
   int nAs = 0;
   int nCs = 0;
   int nGs = 0;
   int nTs = 0;
   int nNs = 0;
   int nOthers = 0;
   int nNOrLowercase = 0;
   long long llNumberOfLowercaseBasesAllSequences = 0;
   int nCsJustOneSequence = 0;
   int nGsJustOneSequence = 0;
   int nNsJustOneSequence = 0;


   if ( argc > 1 ) {
      fprintf( stderr, "Fatal error: reads from stdin\n" );
      return( 1 );
   }


   // needed for printing commas in a number
   setlocale(LC_NUMERIC, "");


   bool bFirstTime = true;
   while( 1 ) {
      c = getchar();
      if ( c == EOF  || c == '>' ) {
         if ( bFirstTime ) {
            bFirstTime = false;
         }
         else {
            szLastHeader[ nHeaderLength ] = 0;
            printf( "%s has %d or %'d bases GC content %.1f%% N's %d\n",  
                    szLastHeader, nNumberOfBases, nNumberOfBases, ( nCsJustOneSequence + nGsJustOneSequence ) * 100.0 / nNumberOfBases, nNsJustOneSequence );
            llNumberOfBasesAllSequences += nNumberOfBases;


            nCsJustOneSequence = 0;
            nGsJustOneSequence = 0;
            nNsJustOneSequence = 0;


            if ( c == EOF ) break;
         }
         bInHeader = true;
         szLastHeader[0] = 0;
         nHeaderLength = 0;
         nNumberOfBases = 0;
      }
      else {
         if ( bInHeader ) {
            if ( c == '\n' ) {
               bInHeader = false;
            }
            else {
               szLastHeader[ nHeaderLength ] = c;
               ++nHeaderLength;
            }
         }
         else {
            if ( !isspace( c )  ) {
               ++nNumberOfBases;
               if ( islower(c) ) {
                  ++llNumberOfLowercaseBasesAllSequences;
               }
               char cUpper = toupper( c );
               if ( cUpper == 'A' ) {
                  ++nAs;
               }
               else if ( cUpper == 'C' ) {
                  ++nCs;
                  ++nCsJustOneSequence;
               }
               else if ( cUpper == 'G' ) {
                  ++nGs;
                  ++nGsJustOneSequence;
               }
               else if ( cUpper == 'T' ) {
                  ++nTs;
               }
               else if ( cUpper == 'N' ) {
                  ++nNs;
                  ++nNsJustOneSequence;
               }
               else {
                  ++nOthers;
               }

               if ( islower(c) || ( c == 'N' ) ) {
                  ++nNOrLowercase;
               }
            }
         }
      }
   }

   printf( "number of bases all sequences: %lld %'lld\n", 
        llNumberOfBasesAllSequences, 
        llNumberOfBasesAllSequences );
   printf( "lowercase bases all sequences: %lld %'lld ( %.2f %% )\n",
           llNumberOfLowercaseBasesAllSequences,
           llNumberOfLowercaseBasesAllSequences,
           llNumberOfLowercaseBasesAllSequences * 100.0 / llNumberOfBasesAllSequences );

   printf( "A: %'d C: %'d G: %'d T: %'d N: %'d (%.1f %%) Other: %'d\n", nAs, nCs, nGs, nTs, nNs, nNs * 100.0 / llNumberOfBasesAllSequences , nOthers );

   printf( "GC content: %.1f%%\n", ( nCs + nGs ) * 100.0 / llNumberOfBasesAllSequences );

   printf( "N or lowercase: %d %'d\n",  nNOrLowercase,  nNOrLowercase );

   return 0;
}


