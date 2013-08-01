#include  <stdio.h>
#include <string.h> 
     #include  <math.h>
     #include  "fold_vars.h"
     #include  "fold.h"
 #include  "utils.h"

 void main(int argc, char *argv[] )
     {
        char *seq1=argv[1], *struct1;
        float e;

        temperature = 37.;
        
        do_backtrack = 0;    
        update_fold_params();
        //initialize_fold(strlen(seq1));
        

        if(argc==2){ 
          struct1 = (char* ) space(sizeof(char)*(strlen(seq1)+1));  
          e =  fold(seq1, struct1);
        }
        else if(argc==3){
          dangles = 0;
          update_fold_params();
          struct1 = argv[2];
          e = energy_of_struct(seq1,struct1);
        }
        else{
          printf("usage: dG <sequence> [structure]\n");
          exit(0);
        }
     
        printf("%f",e);
     
        free_arrays();     /* free arrays used in fold() */
}

