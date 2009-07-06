/*  <Multikulti algorithm: an island model implementation of a GA> */
/*     Copyright (C) <2009>  <Lourdes Araujo, Juan Julian Merelo> */

/*     This program is free software: you can redistribute it and/or modify */
/*     it under the terms of the GNU General Public License as published by */
/*     the Free Software Foundation, either version 3 of the License, or */
/*     (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <time.h>  
#include <iomanip> // for setprecision
//#include "/usr/share/pvm3/include/pvm3.h"
#include "pvm3.h"
#include "funciones.h"

using namespace std;

//#define FUNCIONES "/home/lurdes/sistemas/multikulty_var/funciones"
#define FUNCIONES "/home/juanlu/lurdes/multikulti_pvm/funciones"

#define LON_CAB    3
#define N_WAM      3
#define LON_MAX_MEN 500
#define N_MAX_MEN   12

/* datos externos */
//extern struct TpRegla Reglas[]; 
void Visualizar_chromosome(TIndividuo indi);

/*{{{  static procedures */
main(int argc, char *argv[] )   
{
    int mytid;                  /* my task id */
	int tids[32];				/* slave task ids */
    int n, nproc, numt, i, who, msgtype, nhost, narch;
    float data[100], result[32];
    struct pvmhostinfo *hostp[32];
    int longi;
    char mensaje[LON_MAX_MEN];
    int j,k,l,dir;
    char gen;
    /*---------------------------------------------------------------------*/
    /* estructura de datos para los mejores cromosomas de los analizadores */
    /*---------------------------------------------------------------------*/
    TPoblacion mejor_chro; 
    int *num_ites;  
    int mejor; // posicion del mjor recibido
    float fitness_mejor; // fitness del mejor recibido
    char bit_char; // para copia menesajes
    /*-----------------------------------*/
    /* datos para las medidas de tiempos */
    /*-----------------------------------*/
    struct timeval tv1, tv2;
    int tiempo;   /* en us */     

    /* enroll in pvm */
    mytid = pvm_mytid();

    /* Para que las tareas hijas escriban en salida standard */
    pvm_catchout(stdout);

    /*--------------------------*/
    /* creacion de tareas hijas */
    /*--------------------------*/
    /* Set number of slaves to start */
    /* Can not do stdin from spawned task */
    //puts("Cuantos programas esclavos (1-32)?");
    //scanf("%d", &nproc);
//TEMPORAL
nproc = 8;

    /*--------------------------------------------------*/
    /* se reserva espacio para los cromosomas a recibir */
    /*--------------------------------------------------*/
/*
    if ((mejor_chro = new TIndividuo[nproc])==NULL)
    {
      printf("ERROR en new de mejor_chro\n");
      exit(1);
    }     
*/
    /*---------------------------------------------------------*/
    /* se reserva espacio para el num de iteraciones a recibir */
    /*---------------------------------------------------------*/
    if ((num_ites = new int[nproc])==NULL)
    {
      printf("ERROR en new de num_ites\n");
      exit(1);
    }        
    /*-------------------------------------------*/
    /* se envia mensaje de inicio a tareas hijas */
    /*-------------------------------------------*/
/*
    pvm_initsend(PvmDataDefault);

    pvm_pkint(&nproc, 1, 1);
    pvm_pkint(tids, nproc, 1);
    pvm_mcast(tids, nproc, 0);
*/

    gettimeofday(&tv1, (struct timezone*)0);  
    /* start up slave tasks */
    //numt=pvm_spawn(FUNCIONES, (char**)0, 0, "", nproc, tids);
    numt=pvm_spawn(FUNCIONES, argv+1, 0, "", nproc, tids);
//cout << "numt = " << numt << endl;
/*
    numt=pvm_spawn(FUNCIONES, (char**)0, 1, "bach", nproc, tids);
*/

    if( numt < nproc ){
       printf("Trouble spawning slaves. Aborting. Error codes are:\n");
       for( i=numt ; i<nproc ; i++ ) {
          printf("TID %d %d\n",i,tids[i]);
       }
       for( i=0 ; i<numt ; i++ ){
          pvm_kill( tids[i] );
       }
       pvm_exit();
       exit(1);
    }
    
    /*------------------------------*/
    /* Wait for results from slaves */
    /*------------------------------*/
    msgtype = 5;
//cout << "CONTROL: esperando mensajes " << endl;
    mejor = 0;
    fitness_mejor = 0;
    for( i=0 ; i<nproc ; i++ )
    {
       pvm_recv( -1, msgtype );
       pvm_upkint( &who, 1, 1 );
       pvm_upkint(num_ites+who, 1, 1 );   
       pvm_upkint(&longi, 1, 1 );
//printf("CONTROL: recibido mensaje de FUNCIONES %d longi = %d \n",who, longi);
//cout << "NUM_VAR = " << NUM_VAR << endl;
//cout << "longi = " << longi << endl;
       /*-------------------------------------------------*/
       /* se reserva espacio para os individuos a recibir */
       /*-------------------------------------------------*/
       TIndividuo indiv;
       TGenes genes;
       for(l=0; l < nproc; l++){
         indiv.vec_genes.clear();
         for(k=0; k< NUM_VAR; k++){
           genes.clear();
           //for(j=0; j < longi; j++){
           for(j=0; j < LON_CAD; j++){
             gen = 0;
             genes.push_back(gen);
           }
           indiv.vec_genes.push_back(genes);
         }
         mejor_chro.push_back(indiv);
       }
       //pvm_upkbyte((char *)(mejor_chro+who),longi,1);
       // se reciben los genes
       for(k=0; k < NUM_VAR; k++){
         for(j=0; j < longi; j++){
           pvm_upkbyte(&bit_char,1,1);
//cout << "CONTROL: recibido mensaje de FUNCIONES " << who << "bit = " << bit_char<< endl;
/*
           if (bit_char == '0')
             gen = 0;
           else
             gen = 1;
*/
           //mejor_chro[who].vec_genes[k].push_back(gen);
           mejor_chro[who].vec_genes[k].push_back((bool)bit_char);
         }
       }

//printf("CONTROL: recibido mensaje de FUNCIONES %d \n",who);
       // se recibe la aptitud
       pvm_upkdouble(&(mejor_chro[who].aptitud), 1, 1 );
//---------------//
// MMDP y P-peaks//
//---------------//
       /* se comprueba cual es el mejor de los recibidos */
//cout << "CONTROL: fitness  de " << who << " = " <<  mejor_chro[who].aptitud << endl;
//cout << "CONTROL: fitness_mejor = " << fitness_mejor << endl;
       if (mejor_chro[who].aptitud > fitness_mejor){
         mejor = who;
         fitness_mejor = mejor_chro[who].aptitud;
       }
       /* se imprime el mensaje */
/*
printf("!!!!!!!!!!!!!!EL MEJOR de %d ES:\n ", who);
    Visualizar_chromosome(mejor_chro[who]);
*/
//printf("NUM ITERACIONES = %d\n", nite);
//printf("DESVIACION = %f\n", desviacion);   
// TERMINACION
if (mejor_chro[who].aptitud == NUM_VAR)
  break;
    }
    /*----------------------*/
    /* se calcula el tiempo */
    /*----------------------*/
    gettimeofday(&tv2, (struct timezone*)0);
    tiempo = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
                                                                    
    /* se imprime el mejor */
//cout << "!!!!!!!!!!!!!!EL MEJOR ES:\n " << endl;
    //Visualizar_chromosome(mejor_chro[mejor]);
cout << "NUM_ITES = " << num_ites[mejor] << endl;
cout << "TIEMPO (en micross.) = " << tiempo << endl;
          

    /* Program Finished exit PVM before stopping */
/*
       for( i=0 ; i<nproc ; i++ )
          pvm_kill( tids[i] );
*/
    /*-------------------------------------------------------*/
    /* se envia el mensaje de terminacion a las tareas hijas */
    /*-------------------------------------------------------*/
    for( i=0 ; i<nproc ; i++ ){
      pvm_initsend(PvmDataDefault);
      msgtype = 6; // MEN_FIN
      pvm_send( tids[i], msgtype );
    }


    /*-------------------------------------------------*/
    /* se recibe el numero de evaluaciones de cada una */
    /*-------------------------------------------------*/
    long num_eva_total = 0;
    int num_eva;

    msgtype = 7; // MEN_NUM_EVA
    for( i=0 ; i<nproc ; i++ ){
      // se recibe el numero de evaluaciones 
       pvm_recv( -1, msgtype );
       pvm_upkint( &num_eva, 1, 1 );
//cout << "num_eva de " << i << " = " << num_eva << endl;
       num_eva_total = num_eva_total + num_eva;
    }
    cout << "num_eva_total = " << num_eva_total << endl;

       for( i=0 ; i<nproc ; i++ )
          pvm_kill( tids[i] );

    pvm_exit();
    exit(1);
}


void Visualizar_chromosome(TIndividuo indi)
{
  int i,j;

  cout << "INDIVIDUO:" << endl;
  cout << setprecision(10); // show 10 digits 
  cout << "Aptitud = " << indi.aptitud << endl;
  for(j=0; j < NUM_VAR; j++){
    for (i=0;i<indi.vec_genes.size(); i++)
      cout << indi.vec_genes[j][i] << " ";
    cout << endl;
  }
}


