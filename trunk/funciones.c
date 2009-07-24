 <Multikulti algorithm: an island model implementation of a GA>
    Copyright (C) <2009>  <Lourdes Araujo, Juan Julian Merelo>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fstream>
#include "funciones.h"

//int NUM_VAR = 20; // se pasa funciones.h para que la vea control
Tvec_genes picos;  // se pone como global para que la acceda aptitud
int num_eva = 0;
//TCache cache;

int main( int argc, char *argv[] )    
{
  int tam_pob;           // tamanho de poblacion
  int lcrom;             // longitud de las cadenas
  int generacion;        // contador de generaciones
  int num_gen;           // numero de generaciones
  TPoblacion pob;
  int semilla;
  int pos_mejor;
  float sumaptitud;
  float prob_cruce, prob_mut;
  //TVector evo;  // vector con datos de evolucion
  //TDato dato;
  Tvec_int orden_pob;    // vector con las posicones de la pob ordenadas por aptitud

  //double x_min = -40.0;
  //double x_max = 40.0;
  //double x_min = 0.0;
  //double x_max = 10.0;
  double x_min = 0.0;
  double x_max = 1.0;

  const float TOL = 0.0000001;
  /*---------------------------*/
  /* Datos para comunicaciones */
  /*---------------------------*/
  int mytid;       /* my task id */
  int *tids;    /* task ids   */
      int n, me, i, j, k, nproc, master, msgtype;
      float data[100], result;
  int longi;
  int who;
  int per_inter; // porcentanje de poblacion que se intercambia
  int intervalo; // intervao entre intercambios
  int num_inter; // numero de intercambios
  int seleccionado; // individuo seleccionado para algun intercambio
  int receptor; // tid del proceso al que se envia el mensaje
  int* sustituidos; // array de individuos seleccionados para su sustitucion
  int num_sustituidos; // contador de individuos sustituidos
  int sustituido;
  bool encontrado; // auxiliar para la busqueda de sustituidos
  char bit_char;
  char** representa_genotipo;
  char** secuencia_consenso;
  int reduccion_elite;      // pala la seleccion de indis a enviar

  /*------------------------------------*/
  /* inicilizacion de la comunicaciones */
  /*------------------------------------*/
  /* enroll in pvm */
      mytid = pvm_mytid();

  //lcrom = (int)ceil(log(1 + (x_max - x_min)/TOL) / log((double)2));
  //lcrom = LON_CAD;
//cout << "lcrom = " << lcrom << endl;

  /*-----------------------------*/
  /* inicializaci'on de la cache */
  /*-----------------------------*/
  //cache.clear();

  //--------------------------//
  // P-Peaks P=100, N=LON_CAD //
  //--------------------------//
/*
  picos.clear();
  // se generan 100 cadenas aleatorias de 64 bits.
  srand(50);  // se usa la misma semilla en todos los procesadores
  //Tvec_genes picos;  // se pone como global para que la acceda aptitud
  TGenes genes_picos;
  int gen;
  float f;
  for(i=0; i < NUM_PICOS; i++){
    genes_picos.clear();
    for(j=0; j < LON_CAD; j++){
      f = rand();
      if (f < RAND_MAX / 2 )
        gen = 0;
      else
        gen = 1;
      genes_picos.push_back(gen);
    }
    picos.push_back(genes_picos);
  }
  lcrom = LON_CAD;
//cout << "MAIN: generados picos " << endl;
*/
  //--------------------------//
  // MMDP                     //
  //--------------------------//
  lcrom = LON_CAD;

  /*--------------------------------------*/
  /* preparacion de datos para multikulti */
  /*--------------------------------------*/
  representa_genotipo = new char*[NUM_VAR];
  for(i=0; i < NUM_VAR; i++)
    representa_genotipo[i] = new char[lcrom];
  secuencia_consenso = new char*[NUM_VAR];
  for(i=0; i < NUM_VAR; i++)
    //secuencia_consenso[i] = new char[lcrom];
    secuencia_consenso[i] = new char[LON_CAD];
  
  //num_inter = 1;
  //sustituidos = new int[num_inter];
  /*-----------------------------------------------*/
  /* Preparacion para la funcion de prueba P-peaks */
  /*-----------------------------------------------*/


  // inicializaciones
  //srand(semilla);
  time_t t;
  srand((unsigned) time(&t));


  // TAMAÑOS PEQUES
  //tam_pob = 256;           // tamanho de poblacion (1 proc)
  //tam_pob = 128;           // tamanho de poblacion (2 procs)
  //tam_pob = 64;           // tamanho de poblacion (4 procs) 
  //tam_pob = 32;           // tamanho de poblacion (8 procs)

  // TAMANO GRANDES
  //tam_pob = 512;           // tamanho de poblacion (1 proc)
  //tam_pob = 256;           // tamanho de poblacion (2 procs)
  //tam_pob = 128;           // tamanho de poblacion (4 procs) 
  tam_pob = 64;           // tamanho de poblacion (8 procs)

  // DEMASIADO GRANDES
  //tam_pob = 2048;           // tamanho de poblacion (1 proc)
  //tam_pob = 1024;           // tamanho de poblacion (2 procs)
  //tam_pob = 512;           // tamanho de poblacion (4 procs) 
  //tam_pob = 256;           // tamanho de poblacion (8 procs)

  //num_gen = 100000;          // numero de generaciones
  num_gen = 3000;          // numero de generaciones
  //prob_cruce= (float)20 / (float)100;
  //prob_mut = 10 / (float)100;

  prob_cruce= (float)40 / (float)100;  // P-Peaks y mmdp
  prob_mut = 1 / (float)100;  // P-Peaks y mmdp

  //prob_cruce= (float)40 / (float)100; //wP-Peaks
  //prob_mut = 10 / (float)100; //wP-Peaks

  //prob_cruce= (float)60 / (float)100;  // va muy bien
  //prob_mut = 10 / (float)100;  //va muy bien

  per_inter = 20;
  intervalo = 20;
  orden_pob.clear();
  for(i=0; i < tam_pob; i++)
    orden_pob.push_back(i);

//PVM
  /*--------------------------*/
  /* Receive data from master */
  /*--------------------------*/
/*
  msgtype = 0;
  pvm_recv( -1, msgtype );
      pvm_upkint(&nproc, 1, 1);
      pvm_upkint(tids, nproc, 1);
*/
       /* determine the size of my sibling list */

  tids = new int[20];
        nproc = pvm_siblings(&tids);

  /* Determine which slave I am (0 -- nproc-1) */
  for( i=0; i<nproc ; i++ )
     if( mytid == tids[i] ){ me = i; break; }

//cout << "FUNCIONES " << me << " recibido mensaje de CONTROL" << endl;


  /*---------------------------------------------*/
  /* preparacion del nombre del archivo entropia */
  /*---------------------------------------------*/
  char sme[4];
  string nomb_arch_entropia;
  nomb_arch_entropia  = "/home/juanlu/lurdes/multikulti_pvm/datos/entropia_n";
  sprintf(sme,"%d",me);
  nomb_arch_entropia.append(sme);
  nomb_arch_entropia.append("_mk_cons_elite8.dat");
  //nomb_arch_entropia.append("_best.dat");
cout << " nomb_arch_entropia  = " <<  nomb_arch_entropia << endl;
  /*-------------------------------------*/
  /* se prepara el archivo de entropia */
  /*-------------------------------------*/
  ofstream fentropia(nomb_arch_entropia.c_str());
  if (!fentropia){
    cout << "ERROR en creacion de archivo de entropia " << nomb_arch_entropia << endl;
    exit(1);
  }

  /*-----------------------------------------*/
  /* semilla personalizada para cada proceso */
  /*-----------------------------------------*/
  /* srandom(1); */
  /* srandom(times()); */
  //sscanf(argv[8], "%d", &semilla );
  //semilla = 55;
  semilla = (unsigned) time(&t);
  if (me > 0){
    //semilla = semilla * (me+1);
    //semilla = (long)pow((double)semilla, (double)me+1);
    semilla = (int)(semilla/(double)(me+1));
//cout << "semilla  = " << semilla << endl;
    //semilla = me;
  }

  srandom(semilla); // para ejecutar en varios procesadores
  srand(semilla); // para ejecutar en varios procesadores

  // se genera la poblacion inicial
  pob = poblacion_inicial(tam_pob, lcrom, x_min, x_max);
//cout << "FUNCIONES " << me << " creada POB INI" << endl;
  // se calculan las estadisticas
  evaluacion(pob, tam_pob, pos_mejor, sumaptitud, orden_pob);
//cout << "FUNCIONES " << me << " evaluada POB INI" << endl;

  // se entra en el bucle de evolucion
  for (generacion=0; generacion < num_gen; generacion++){
    // se seleccionan tam_pob individuos con oportunidades de reproducirse
    seleccion(pob, tam_pob);

    reproduccion(pob, tam_pob, lcrom, prob_cruce, x_min, x_max);

    // se aplica el operador de mutacion
    mutacion(pob, tam_pob, lcrom, prob_mut, x_min, x_max);

    evaluacion(pob, tam_pob, pos_mejor, sumaptitud, orden_pob);
//if (generacion % 1000 == 0){
//cout << "GENERACION " << generacion << endl;
//cout << "Aptitud del mejor(" << pos_mejor << ") = " << pob[pos_mejor].aptitud << " en " << pob[pos_mejor].vec_x[0] << endl;
//cout << "Aptitud media = " << sumaptitud/tam_pob << endl;
//}

//cout << "me = " << me << " entropia = " << calcula_entropia(pob, tam_pob, lcrom) << endl;
    if (generacion % (intervalo+1) == 0)
      fentropia << calcula_entropia(pob, tam_pob, lcrom) << endl;

    /*----------------------------------------------------------------------*/
    /* se comprueba si se ha recibido algun mensaje de terminacion (tipo 6) */
    /*----------------------------------------------------------------------*/
    msgtype = 6;
    if (pvm_probe(-1,msgtype) != 0){ // MEN_TERMINA
break;
//cout << " PROC " << me << " RECIBE mensaje MEN_TERMINA " << endl;
      /*------------------------------------------------------------------------*/
      /* se recibe un genotipo que representa a la poblacion del proc. anterior */
      /*------------------------------------------------------------------------*/
      msgtype = 6;
      pvm_recv( -1, msgtype );

      /*----------------------------------------------*/
      /* se envia el numero de evaluaciones a control */
      /*----------------------------------------------*/
      pvm_initsend( PvmDataDefault );
      // se envia el numero de evaluaciones
      pvm_pkint( &num_eva, 1, 1 );

      msgtype = 7; // MEN_NUM_EVA
      master = pvm_parent();
      pvm_send( master, msgtype );
      //pvm_exit();
      //exit(1);
    }


    /*------------------------------------------------*/
    /* envio y recepcion de individuos de otros demes */
    /*------------------------------------------------*/
    if ( (nproc > 1) && (generacion % intervalo == 0) ){
    //if (generacion % intervalo == 0) {
//cout << " PROC " << me << " entra en el IF de mensajes "  << endl;
      /*------------------------*/
      /* se inicializa el envio */
      /*------------------------*/
      pvm_initsend( PvmDataDefault );
      /*-------------------------------------------------------------*/
      /* se envia una representacion del genotipo de la poblacion al */
      /* procesador anterior                                         */
      /*-------------------------------------------------------------*/
      /*-------------------------------------*/
      /* se copia el cromosoma en el mensaje */
      /*-------------------------------------*/
      calcula_secuencia_consenso(pob, tam_pob, secuencia_consenso, lcrom);    // la secuencia consenso
      for(i=0; i < NUM_VAR; i++){
        for (j =0; j< LON_CAD; j++){
          bit_char = secuencia_consenso[i][j];    // la secuencia consenso
          pvm_pkbyte(&bit_char, 1, 1 );
        }
      }
      /*---------------------------*/
      /* se envia el mejor         */
      /*---------------------------*/
/*
      for(i=0; i < NUM_VAR; i++){
        for (j =0; j< LON_CAD; j++){
          bit_char = pob[pos_mejor].vec_genes[i][j];  // el mejor
          pvm_pkbyte(&bit_char, 1, 1 );
        }
      }
*/
      /*---------------------------*/
      /* se envia proceso anterior */
      /*---------------------------*/
      msgtype = 4; // MEN_GENOTIPO
      //receptor = (me + 1) % nproc;
      if (me - 1 < 0)
        receptor = nproc - 1;
      else
        receptor = me - 1;
//cout << " PROC " << me << " envia mensaje MEN_GENOTIPO a PROC " << receptor << endl;
      pvm_send( tids[receptor], msgtype );

      msgtype = 4;
      if (pvm_probe(-1,msgtype) != 0){ // MEN_GENOTIPO
//cout << " PROC " << me << " RECIBE mensaje MEN_GENOTIPO de PROC " << receptor << endl;
        /*------------------------------------------------------------------------*/
        /* se recibe un genotipo que representa a la poblacion del proc. anterior */
        /*------------------------------------------------------------------------*/
        msgtype = 4;
        pvm_recv( -1, msgtype );
        // se copian los genes
        for(i=0; i < NUM_VAR; i++){
          for (j =0; j< LON_CAD; j++){
            pvm_upkbyte(&bit_char,1,1);
            representa_genotipo[i][j] = bit_char;
          }
        }
        /*-----------------------------------------------------*/
        /* se envia el individuo mas diferente de la poblacion */
        /*-----------------------------------------------------*/

        // se envia al siguiente proc
        /*------------------------*/
        /* se inicializa el envio */
        /*------------------------*/
        pvm_initsend( PvmDataDefault );
        /*----------------------------------------------------------*/
        /* se envian tam_pob * per_inter individuos seleccionados */
        /* proporcionalmente a su fitness                           */
        /*----------------------------------------------------------*/
        //num_inter = tam_pob * per_inter / 100;
        num_inter = 1;
        longi = num_inter * sizeof(lcrom);
        pvm_pkint( &longi, 1, 1 );
  
        longi = 0;
        for (i=0; i< num_inter; i++){
          //seleccionado = pos_mejor;
          //seleccionado = seleccion_emigrante_mejor(pob, tam_pob, sumaptitud); // con sesgo
          //seleccionado = alea(0, tam_pob-1);
          //seleccionado  = busca_diferente(pob, tam_pob, representa_genotipo, lcrom);

          //reduccion_elite = 1;
          //reduccion_elite = 2;
          //reduccion_elite = 4;
          //reduccion_elite = 8;
          //reduccion_elite = 16;
          //reduccion_elite = 32;
          //seleccionado  = busca_diferente_elite(pob, tam_pob, representa_genotipo, lcrom, orden_pob, reduccion_elite);

          // con elite fija
          //seleccionado  = busca_diferente_elite_fija(pob, tam_pob, representa_genotipo, lcrom, orden_pob, 4);
          seleccionado  = busca_diferente_elite_fija(pob, tam_pob, representa_genotipo, lcrom, orden_pob, 8);
          //seleccionado  = busca_diferente_elite_fija(pob, tam_pob, representa_genotipo, lcrom, orden_pob, 16);
          /*-------------------------------------*/
          /* se copia el cromosoma en el mensaje */
          /*-------------------------------------*/
          for(k=0; k < NUM_VAR; k++){
            for (j =0; j< lcrom; j++){
              bit_char = pob[seleccionado].vec_genes[k][j];
              pvm_pkbyte(&bit_char, 1, 1 );
//if( bit_char == 1)
//cout << '1';
//else if (bit_char == 0)
//cout << '0';
            }
          }
// se envia la aptitud
pvm_pkdouble( &(pob[seleccionado].aptitud), 1, 1 );
          //longi = longi + NUM_VAR*lcrom;
          //longi = longi + NUM_VAR*LON_CAD;
          longi = longi + NUM_VAR*LON_CAD + sizeof(double);
        }
        /*-------------------------------*/
        /* se envia al siguiente proceso */
        /*-------------------------------*/
        msgtype = 3; // MEN_INTER
        receptor = (me + 1) % nproc;
//cout << " PROC " << me << " envia mensaje MEN_INTER a PROC " << receptor << endl;
//cout << " PROC " << me << " fitness de enviado = " << pob[seleccionado].aptitud << endl;
        pvm_send( tids[receptor], msgtype );
      }
      if (pvm_probe(-1,3) != 0){ // MEN_INTER
        /*------------------------------------------------------------*/
        /* se reciben otros tantos individuos, que sustituyen a otros */
        /* seleccionados aleatoriamente                               */
        /*------------------------------------------------------------*/
        msgtype = 3;
        pvm_recv( -1, msgtype );
        pvm_upkint(&longi, 1, 1);
//cout << " PROC " << me << " recibido mensaje MEN_INTER" << endl;
        /*-------------------------------------------------------------------------*/
        /* los cromosomas recibidos se copian e otros seleccionados aleatoriamente */
        /*-------------------------------------------------------------------------*/
/*
        num_sustituidos = 0;
        //sustituidos = new int[num_inter];
        while (num_sustituidos < num_inter){
          sustituido = alea(0, tam_pob-1);
          // se comprueba que no estuviera ya seleccioando 
          j = 0;
          encontrado = false;
          while( (j < num_sustituidos) && !encontrado){
            if (sustituidos[j] != sustituido)
              j++;
            else
              encontrado = true;
          }
          if (!encontrado){
            sustituidos[num_sustituidos] = sustituido;
            num_sustituidos++;
          }
        }
*/
        /*--------------------------------------------------------------------*/
        /* los cromosomas recibidos se copian en las posiciones de los peores */
        /*--------------------------------------------------------------------*/
        num_sustituidos = 0;
        num_inter=1;
        sustituidos = new int[num_inter];
        while (num_sustituidos < num_inter){
          sustituidos[num_sustituidos] = orden_pob[tam_pob - 1 - num_sustituidos];
          num_sustituidos++;
        }
//cout << "sustituidos[0] = " << sustituidos[0] << endl;
        /*----------------------------*/
        /* se hacen las sustituciones */
        /*----------------------------*/
        for (i =0; i< num_inter; i++){
//cout << " PROC " << me << " se sustituye en indivuo " << sustituidos[i] << " de aptitud " << pob[sustituidos[i]].aptitud << endl;
          // se copian los genes
          for(k=0; k < NUM_VAR; k++){
            for (j =0; j< LON_CAD; j++){
              pvm_upkbyte(&bit_char,1,1);
              pob[sustituidos[i]].vec_genes[k][j] = bit_char;
/*
if( bit_char == 1)
cout << '1';
else if (bit_char == 0)
cout << '0';
*/
//cout << pob[sustituidos[i]].vec_genes[k][j]; 
            }
          }
//cout << endl;
          //pob[sustituidos[i]].aptitud = adaptacion(pob[sustituidos[i]], lcrom, x_min, x_max);
pvm_upkdouble(&(pob[sustituidos[i]].aptitud), 1, 1 );
//cout << " PROC " << me << " aptitud del sustituto =  " << pob[sustituidos[i]].aptitud << endl;
//cout << " PROC " << me << " aptitud del sustituto =  " << adaptacion(pob[sustituidos[i]], lcrom, x_min, x_max) << endl;
        }

        /*-------------------------*/
        /* se liberan los recursos */
        /*-------------------------*/
        delete [] sustituidos;
      } // pvm_probe
    }
//cout << " PROC " << me << " ha pasado el IF de mensajes "  << endl;
if(pob[pos_mejor].aptitud == NUM_VAR)  //MMDP y P-peaks
break;
  }
  /*---------------------------------------*/
  /* se envia el mejor cromosoma a control */
  /*---------------------------------------*/
//cout << "num_eva[" << me << "] = " << num_eva << endl;
//cout << "PROC " << me << " Se va a enviar mejor cro a CONTROL " << endl;
  pvm_initsend( PvmDataDefault );
  pvm_pkint( &me, 1, 1 );
  pvm_pkint( &generacion, 1, 1 );
  //longi = sizeof(TIndividuo);
  longi = lcrom;
  pvm_pkint( &longi, 1, 1 );
  //pvm_pkbyte(mensaje, longi, 1 );
  //pvm_pkbyte((char *)(pob+pos_mejor), sizeof(TIndividuo), 1 );
//cout << "PROC " << me << " Se va a enviar los genes " << endl;
  for(i=0; i < NUM_VAR; i++){
    for (j =0; j< lcrom; j++){
      bit_char = pob[pos_mejor].vec_genes[i][j];
/*
      if (pob[pos_mejor].vec_genes[i][j] == true)
        bit_char = '1';
      else if (pob[pos_mejor].vec_genes[i][j] == false)
        bit_char = '0';
      else
        bit_char = '2';
*/
//cout << "PROC " << me << " Se va a enviar el gen " << bit_char << endl;
      pvm_pkbyte(&bit_char, 1, 1 );
    }
  }

  // se envia la adaptacion
  pvm_pkdouble( &(pob[pos_mejor].aptitud), 1, 1 );

  msgtype = 5; // MEN_CHRO
  master = pvm_parent();
  pvm_send( master, msgtype );
//cout << "PROC " << me << " enviado mejor cro a CONTROL " << endl;

  /*---------------------*/
  /* se liberan recursos */
  /*---------------------*/
  delete [] representa_genotipo;

  /*---------------------------------------------------------*/
  /* se espera a que Control envia la se\~nal de terminacion */
  /*---------------------------------------------------------*/
  while(1){
    /*----------------------------------------------------------------------*/
    /* se comprueba si se ha recibido algun mensaje de terminacion (tipo 6) */
    /*----------------------------------------------------------------------*/
    msgtype = 6;
    if (pvm_probe(-1,msgtype) != 0){ // MEN_FIN
//cout << "Proc " << me << " recibido MEN_FIN " << endl;
//cout << "se envia a CONTROL num_eva " << me << num_eva << endl;
      /*------------------------------------------------------------------------*/
      /* se recibe un genotipo que representa a la poblacion del proc. anterior */
      /*------------------------------------------------------------------------*/
      msgtype = 6;
      pvm_recv( -1, msgtype );

      /*----------------------------------------------*/
      /* se envia el numero de evaluaciones a control */
      /*----------------------------------------------*/
      pvm_initsend( PvmDataDefault );
      // se envia el numero de evaluaciones
      pvm_pkint( &num_eva, 1, 1 );

      msgtype = 7; // MEN_NUM_EVA
      master = pvm_parent();
      pvm_send( master, msgtype );
      //pvm_exit();
      //exit(1);
    }
  }
  while(1)
    pvm_recv( -1, -1 );

}
/************************************************************/
/*                                                          */
/* funcion: poblacion_inicial                               */
/* Genera la poblacion inicial aleatoriamente               */
/*                                                          */
/************************************************************/
TPoblacion poblacion_inicial(int tam_pob, int lcrom, float x_min, float x_max)
{
  TPoblacion pob;
  TIndividuo indiv;
  int i;
  for(i=0; i < tam_pob; i++){
    indiv = genera_indiv(lcrom, x_min, x_max);
    pob.push_back(indiv);
  }
  return pob;
}
/************************************************************/
/*                                                          */
/* funcion: genera_indiv                                    */
/* Genera un individuo de la poblacion inicial              */
/*                                                          */
/************************************************************/
TIndividuo genera_indiv(int lcrom, float x_min, float x_max)
{
  int i,j;
  TIndividuo indiv;
  TGenes genes;
  int gen;
  float f;

  indiv.vec_genes.clear();
  /*---------------------------------------------*/
  /* se generan los genotipos para cada variable */
  /*---------------------------------------------*/
  for(i=0; i< NUM_VAR; i++){
    genes.clear();
    for(j=0; j < lcrom; j++){
      f = rand();
      if (f < RAND_MAX / 2 )
        gen = 0;
      else
        gen = 1;
      genes.push_back(gen);
    }
    indiv.vec_genes.push_back(genes);
  }
  indiv.aptitud = adaptacion(indiv, lcrom, x_min, x_max);
  return indiv;
}

/************************************************************/
/*                                                          */
/* funcion: seleccion                                       */
/* Selecciona un unico individuo por el metodo de la ruleta */
/*                                                          */
/************************************************************/
void seleccion(TPoblacion& pob, int tam_pob)
{
  int* sel_super;                  // seleccionados para sobrevivir       
  float prob,f;                    // probabilidad de seleccion
  int pos_super;                   // posicion del superviviente
  TPoblacion pob_aux;              // poblacion auxiliar   
  int i,j;
  //int TAM_ELITE = pob.size()* 2 /100;
  int TAM_ELITE = 1;
  int* sel_elite;

  sel_super = new int[tam_pob];
  sel_elite = new int[TAM_ELITE];
  //------------------------//
  // se selecciona la elite
  //------------------------//
  for (i=0; i < TAM_ELITE; i++){
    sel_elite[i] = 0;
  }
  for (i=0; i < tam_pob; i++){
    pob[i].elite = false;
    j = 0;
    while ( (j < TAM_ELITE) && (pob[sel_elite[j]].aptitud >= pob[i].aptitud)){
      j++;
    }
    if (j < TAM_ELITE){
     sel_elite[j] = i;
    }
  }
  for (i=0; i < TAM_ELITE; i++){
    pob[sel_elite[i]].elite = true;
//cout << " seleccion: apt elite " << i << " = " << pob[sel_elite[i]].aptitud << endl;
  }

  // se seleccionan tam_pob individuos para reproducirse
  // se generan numeros aleatorios entre 0 y 1, seleccionan individuos
  // de acuerdo con su puntuacion acumulada
 
  for (i=0; i < tam_pob; i++){
    if (pob[i].elite == true)  //ELITISMO
      sel_super[i] = i;
    else{
      f = rand();
      prob = (f/(RAND_MAX+1.0));
      pos_super = 0;
      while ((prob > pob[pos_super].punt_acu) && (pos_super < tam_pob))
        pos_super++;
      if (pos_super < tam_pob)
        sel_super[i] = pos_super;
      else
        sel_super[i] = pos_super-1;
    }
  }
  // se genera la poblacion intermedia
  for (i=0; i < tam_pob; i++){
    pob_aux.push_back(pob[sel_super[i]]);
  }
  pob.clear();
  for (i=0; i < tam_pob; i++){
    pob.push_back(pob_aux[i]);
  }                            
  delete [] sel_super;
  delete [] sel_elite;
}

/************************************************************/
/*                                                          */
/* funcion: adaptacion                                      */
/* Evalua la calidad de un individuo                        */
/*                                                          */
/************************************************************/
double adaptacion(TIndividuo& individuo, int lcrom, float x_min, float x_max)
{
  int i;
  Tvec_var x; // vector de fenotipos
  double f; // valor de la funcion a optimizar

  //if (buscar_cadena_cache(cache, individuo.vec_genes[0], f))
  //if (buscar_cadena_cache(cache, individuo, f))
  //{
    //return f;
  //}

  x.clear();
  individuo.vec_x.clear();
  //for(i=0; i < NUM_VAR; i++){
    //individuo.vec_x.push_back(decod(individuo.vec_genes[i], lcrom, x_min, x_max));
  //}
  //x = individuo.vec_x;

  //f = 0.5 - ((sin(fabs(x))*sin(fabs(x))-0.5)/(1.0 + 0.001*(x*x*x*x)));

  //f = 20+exp((double)1)-20*exp((double)(-0.2*fabs(x)))-exp((double)(cos(2*M_PI*x)));
  //f = sin(x)/(1+sqrt(x)+cos(x)/(1+x));
  //f = pow(sin(5*M_PI*x), 6);
  //f = (1/(0.5*sqrt(2*M_PI))) * exp(-0.5*pow(((x-5.0)/0.5),2));
  //f = pow( sin(5*M_PI*x) * exp(-2*log((double)2)*pow(((x-0.1)/0.8), 2)),6);
  //f = 1 + cos(x)/(1+ 0.01*x*x);

  //f = x[0] + fabs(sin(32*M_PI*x[0]));
  //-------//
  //  MMDP //
  //-------//
  int cuenta_unos, j;
  f = 0;
  for(i=0; i < NUM_VAR; i++){
    cuenta_unos = 0;
    for(j=0; j < lcrom; j++){
      if (individuo.vec_genes[i][j] == (char)1)
       cuenta_unos++;
    }
    if(cuenta_unos == 0) f = f + 1.0;
    else if(cuenta_unos == 1) f = f + 0.0;
    else if(cuenta_unos == 2) f = f + 0.360384;
    else if(cuenta_unos == 3) f = f + 0.640576;
    else if(cuenta_unos == 4) f = f + 0.360384;
    else if(cuenta_unos == 5) f = f + 0.0;
    else if(cuenta_unos == 6) f = f + 1.0;
  }
  //---------//
  // P-peaks //
  //---------//
/*
//cout << "adaptacion: entre en p-peaks " << endl;
  int min_hamming,hamming,j;
  // se calcula la distancia de Hamming minima a los picos
  min_hamming = LON_CAD;
  for(i=0; i < NUM_PICOS; i++){
    hamming = 0;
    for(j=0; j < LON_CAD; j++){
      if (individuo.vec_genes[0][j] != picos[i][j])
        hamming++;
    }
//cout << "adaptacion: hamming = " << hamming << endl;
    if (hamming < min_hamming)
      min_hamming = hamming;
  }
  f = (float)(LON_CAD - min_hamming)/(float)LON_CAD;
//cout << "adaptacion: sale de p-peaks f = " << f << endl;
*/
  //---------//
  // wP-peaks //
  //---------//
/*
double max_hamming;
double peso;
double hamming,j;
//cout << "adaptacion: entre en p-peaks " << endl;
  // se calcula la distancia de Hamming minima a los picos
  max_hamming = 0;
  for(i=0; i < NUM_PICOS; i++){
    hamming = 0;
    for(j=0; j < LON_CAD; j++){
      if (individuo.vec_genes[0][j] != picos[i][j])
        hamming++;
    }
//cout << "adaptacion: hamming = " << hamming << endl;
    if (i==0) peso=1;
    else peso=0.99;
    hamming = (float)LON_CAD*peso - hamming;
    //if (((double)(LON_CAD*peso) - hamming) > max_hamming)
      //max_hamming = ((double)(LON_CAD*peso) - hamming);
    if (hamming > max_hamming)
      max_hamming =  hamming;
  }
//cout << "adaptacion: max_hamming = " << max_hamming << endl;
  f = (float)max_hamming/(float)LON_CAD;
//cout << "adaptacion: sale de wp-peaks f = " << f << endl;
*/
  
  num_eva++;
  //insertar_cadena_cache(cache, individuo.vec_genes[0], f);
  //insertar_cadena_cache(cache, individuo, f);
  return f;
}
/************************************************************/
/*                                                          */
/* funcion: decod                                           */
/* Decodifica el genotipo de un individuo                   */
/*                                                          */
/************************************************************/
double decod(TGenes genes, int lcrom, float x_min, float x_max)
{
  double x;
  
  x = (float)bin_dec(genes, lcrom) / (pow((double)2,(double)lcrom) - 1);
  x = x_min + (x_max - x_min) * x;
  return x;
}
/************************************************************/
/*                                                          */
/* funcion: bin_dec                                         */
/* Convierte de binario a decimal                           */
/*                                                          */
/************************************************************/
int bin_dec(TGenes genes, int lcrom)
{
  int i, d=0, pot=1;

  for(i=0; i < lcrom; i++){
    //d = d + pot*genes[lcrom - i - 1]; 
    d = d + pot*genes[i]; 
    pot = pot * 2;
  }
  return d;
}

/************************************************************/
/*                                                          */
/* funcion: reproduccion                                    */
/* Selecciona los individuos a reproducirse y aplica el     */
/* operador de cruce                                        */
/*                                                          */
/************************************************************/
void reproduccion(TPoblacion& pob, int tam_pob, int lcrom, float prob_cruce, 
                  float x_min, float x_max)
{
  int sel_cruce[tam_pob];          // seleccionados para reproducirse
  int num_sel_cruce=0;             // num de individuos seleccionados para se cruzados
  float prob;
  int punto_cruce1, punto_cruce2;
  TIndividuo hijo1, hijo2;
  int i;

  // se eligen los individuos a cruzar
  for (i=0; i < tam_pob; i++){
    // se generan tam_pob numeros aleatorios a_i en [0 1)
    prob = (rand()/(RAND_MAX+1.0));   
    // se eligen para el cruce los individuos de las posiciones i con a_i < prob_cruce 
    if (prob < prob_cruce){
      sel_cruce[num_sel_cruce] = i;
      num_sel_cruce++;
    }
  }
  // el numero de selecionados se hace par 
  if ((num_sel_cruce % 2) == 1)
    num_sel_cruce--;

  // se cruzan los individuos elegidos en un punto al azar
  punto_cruce1 = (int)(rand() * (lcrom+1) / (RAND_MAX + 1.0)) + 1;
  for (i=0; i < num_sel_cruce; i=i+2){
    cruce(pob[sel_cruce[i]], pob[sel_cruce[i+1]], hijo1, hijo2, lcrom, punto_cruce1, x_min, x_max);
    //cruce_bipunto(pob[sel_cruce[i]], pob[sel_cruce[i+1]], hijo1, hijo2, lcrom, punto_cruce1, punto_cruce2, x_min, x_max);
    // los nuevos individuos sustituyen a sus progenitores
    //pob[sel_cruce[i]] = hijo1;
    //pob[sel_cruce[i+1]] = hijo2;
    //if ((hijo1.aptitud > adaptacion(pob[sel_cruce[i]], lcrom, x_min, x_max)) || (!pob[sel_cruce[i]].elite)){ //ELITISMO
    if ((hijo1.aptitud > pob[sel_cruce[i]].aptitud) || (!pob[sel_cruce[i]].elite)){ //ELITISMO
      pob[sel_cruce[i]] = hijo1;
    }
    //if ((hijo2.aptitud > adaptacion(pob[sel_cruce[i+1]], lcrom, x_min, x_max)) || (!pob[sel_cruce[i+1]].elite)){
    if ((hijo2.aptitud > pob[sel_cruce[i+1]].aptitud) || (!pob[sel_cruce[i+1]].elite)){
      pob[sel_cruce[i+1]] = hijo2;
    }

  }
}

/************************************************************/
/*                                                          */
/* funcion: cruce                                           */
/* Aplica el operador de cruce                              */
/*                                                          */
/************************************************************/
void cruce(TIndividuo padre1, TIndividuo padre2, TIndividuo& hijo1, TIndividuo& hijo2,
           int lcrom, int punto_cruce, float x_min, float x_max)      
{
  int i,j;
  TGenes genes_hijo1, genes_hijo2;

  hijo1.vec_genes.clear();
  hijo2.vec_genes.clear();

  for(j=0; j < NUM_VAR; j++){
    // primera parte del intercambio: 1 a 1 y 2 a 2
    genes_hijo1.clear();
    genes_hijo2.clear();
    for(i=0; i < punto_cruce; i++){
      genes_hijo1.push_back(padre1.vec_genes[j][i]);
      genes_hijo2.push_back(padre2.vec_genes[j][i]);
    }
    // segunda parte: 1 a 2 y 2 a 1
    for(i=punto_cruce; i < lcrom; i++){
      genes_hijo1.push_back(padre2.vec_genes[j][i]);
      genes_hijo2.push_back(padre1.vec_genes[j][i]);
    }
    hijo1.vec_genes.push_back(genes_hijo1);
    hijo2.vec_genes.push_back(genes_hijo2);
  }
  // se evaluan
  hijo1.aptitud = adaptacion(hijo1, lcrom, x_min, x_max);
  hijo2.aptitud = adaptacion(hijo2, lcrom, x_min, x_max);
}

/************************************************************/
/*                                                          */
/* funcion: cruce                                           */
/* Aplica el operador de cruce                              */
/*                                                          */
/************************************************************/
void cruce_bipunto(TIndividuo padre1, TIndividuo padre2, TIndividuo& hijo1, TIndividuo& hijo2,
           int lcrom, int punto_cruce1, int punto_cruce2, float x_min, float x_max)
{
  int i,j;
  TGenes genes_hijo1, genes_hijo2;

  hijo1.vec_genes.clear();
  hijo2.vec_genes.clear();

  for(j=0; j < NUM_VAR; j++){
    // primera parte del intercambio: 1 a 1 y 2 a 2
    genes_hijo1.clear();
    genes_hijo2.clear();
    for(i=0; i < punto_cruce1; i++){
      genes_hijo1.push_back(padre1.vec_genes[j][i]);
      genes_hijo2.push_back(padre2.vec_genes[j][i]);
    }
    // segunda parte: 1 a 2 y 2 a 1
    for(i=punto_cruce1; i < punto_cruce2; i++){
      genes_hijo1.push_back(padre2.vec_genes[j][i]);
      genes_hijo2.push_back(padre1.vec_genes[j][i]);
    }
    // parte final: 1 a 1 y 2 a 2
    for(i=punto_cruce2; i < lcrom; i++){
      genes_hijo1.push_back(padre1.vec_genes[j][i]);
      genes_hijo2.push_back(padre2.vec_genes[j][i]);
    }
    hijo1.vec_genes.push_back(genes_hijo1);
    hijo2.vec_genes.push_back(genes_hijo2);
  }
  // se evaluan
  hijo1.aptitud = adaptacion(hijo1, lcrom, x_min, x_max);
  hijo2.aptitud = adaptacion(hijo2, lcrom, x_min, x_max);
}

/************************************************************/
/*                                                          */
/* funcion: mutacion                                        */
/* Aplica el operador de mutacion a la poblacion            */
/*                                                          */
/************************************************************/
void mutacion(TPoblacion& pob, int tam_pob, int lcrom, float prob_mut, float x_min, float x_max)    
{
  bool mutado;
  int i,j,k;
  float prob,f;
  double aptitud_old;


  for (i=0; i < tam_pob; i++){
    mutado = false;
    for(k=0; k <NUM_VAR; k++){
      for (j=0; j < lcrom; j++){
        // se genera un numero aleatorio en [0 1)
        f = rand();
        prob = (f/(RAND_MAX+1.0));
        // se mutan aquellos genes con prob < que prob_mut
        //if ((prob < prob_mut) && (!pob[i].elite)){
          //pob[i].vec_genes[k][j] = !(bool)( pob[i].vec_genes[k][j]);
          //mutado = true;
        //}
        if (prob < prob_mut){
          if (!(pob[i].elite)){
            pob[i].vec_genes[k][j] = !(bool)( pob[i].vec_genes[k][j]);
            mutado = true;
          }
/*
          else{  // esta en la elite
            aptitud_old = pob[i].aptitud;
//cout << "MUTACION: aptitud_old = " << aptitud_old << endl;
            pob[i].vec_genes[k][j] = !(bool)( pob[i].vec_genes[k][j]);
            pob[i].aptitud = adaptacion(pob[i], lcrom, x_min, x_max);
            if(aptitud_old > pob[i].aptitud){ // se deshace el cambio
              pob[i].vec_genes[k][j] = !(bool)( pob[i].vec_genes[k][j]);
              pob[i].aptitud = aptitud_old;
            }
            else
              mutado = true;
//cout << "MUTACION: pob[" << i <<"] = " << pob[i].aptitud << endl;
          }
*/
        }
      }
      if (mutado)
        pob[i].aptitud = adaptacion(pob[i], lcrom, x_min, x_max);
    } 
  }
}

/************************************************************/
/*                                                          */
/* funcion: evaluacion                                      */
/* Calcula y actualiza as estadisticas de la poblacion      */
/*                                                          */
/************************************************************/
void evaluacion(TPoblacion& pob, int tam_pob, int& pos_mejor, float& sumaptitud, Tvec_int& orden_pob)
{
  float punt_acu = 0;    // puntacion acumulada de los individuos
  float aptitud_mejor = 0;   // mejor aptitud
  int i,j,k;
  sumaptitud = 0;  // suma de la aptitud

  int pos_orden;

  pos_orden = 0;
  for (i=0; i< tam_pob; i++){
    sumaptitud = sumaptitud + pob[i].aptitud;
    if (pob[i].aptitud > aptitud_mejor){
      pos_mejor = i;
      aptitud_mejor = pob[i].aptitud;
    }
    // se inserta ordenadamente en el vector de orden de la poblacion
    if (i==0)
      orden_pob[0] = 0;
    else{
      j = i-1;
      while((j>=0) && (pob[orden_pob[j]].aptitud < pob[i].aptitud)){
        j--; 
      }
      j++;
      // se desplazan los elementos necesarios para la insercion
      for(k=i-1; k>=j; k--)
        orden_pob[k+1] = orden_pob[k];
      orden_pob[j] = i;
    }
  }
// se visualiza
//for(i=0; i < tam_pob; i++)
//cout << "orden = " << orden_pob[i] << " aptitud = " << pob[orden_pob[i]].aptitud << endl;

  for (i=0; i< tam_pob; i++){
    pob[i].puntuacion = pob[i].aptitud / sumaptitud;
    pob[i].punt_acu = pob[i].puntuacion + punt_acu;
    punt_acu = punt_acu + pob[i].puntuacion;
  }
}

/*******************************************************************/
/*                                                                 */
/* Rutina : Seleccion                                              */
/*                                                                 */
/* Funcion : Selecciona un unico individuo mediante una rueda de   */
/*    la fortuna cuyos sectores son proporcionales a la fitness    */
/*    de los individuos                                            */
/*                                                                 */
/*******************************************************************/
int seleccion_emigrante_mejor(TPoblacion pob, int tam_pob, float sumaptitud)
{
  float punto_alea;  // punto aleatorio de la rueda
  float sum_parcial = 0.0; // suma parcial de fitness de individuos
  int j = 0;

  //punto_alea = (rand() * fabs(sumaptitud)) / (RAND_MAX + 1.0);
  punto_alea = (random() * fabs(sumaptitud)) / (RAND_MAX + 1.0);
//cout << "SELECCION: punto_alea = " << punto_alea << endl;
//char c = getchar();
  do{                // se busca el sector de la rueda que corresponde al punto
    sum_parcial = sum_parcial + pob[j].aptitud;
    j++;
  }
  while ((sum_parcial <= punto_alea) && (j < tam_pob));
//cout << "SELECCION: seleccionado = " << (j - 1) << endl;
//cout << "SELECCION: fitness = " <<  pob[j-1].aptitud << endl;
//char c = getchar();

  return (j-1);
}

/*******************************************************************/
/*                                                                 */
/* Rutina : Alea(min_valor, max_valor)                             */
/*                                                                 */
/* Funcion : Genera un entero aleatoria entre min_valor y max_valor*/
/*                                                                 */
/*******************************************************************/
int alea( int min_valor, int max_valor )
{
  float valor;
/* printf("ALEA: se va a generar un num entre %d y %d\n", min_valor, max_valor);*/
  valor = random();
/* printf("ALEA: se ha genrado %f\n", valor); */

/*
  valor = min_valor + ((int) (max_valor*valor/(RAND_MAX + 1.0)));
*/

  valor = min_valor + ((int) ((max_valor-min_valor+1)*valor/(RAND_MAX + 1.0)));

/* printf("ALEA: se ha convertido a %f\n", valor); */
  return (int)valor;
}

/************************************************************/
/*                                                          */
/* funcion: busca_diferente                               */
/* Busca en la poblaci'on la posicion del individuo mas     */
/* diferente a represeta_genotipo                           */
/*                                                          */
/************************************************************/
int busca_diferente(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom)
{
  int i,j,k;
  int mas_diferente, diferencia, mayor_dif;

  mayor_dif = 0;
  mas_diferente = 0;
  for(i=0; i < tam_pob; i++){
    diferencia = 0;
    for(k=0; k<NUM_VAR;k++){
      for(j=0; j < lcrom; j++){
        if(pob[i].vec_genes[k][j] != *(bool *)(representa_genotipo[k] + j))
          diferencia++;
      }
    }
    if (diferencia > mayor_dif){
      mayor_dif = diferencia;
      mas_diferente = i;
    }
  }
  return mas_diferente;
}
/************************************************************/
/*                                                          */
/* funcion: busca_diferente_medio                           */
/* Busca en la poblaci'on la posicion del individuo con     */
/* una diferencia media del  represeta_genotipo             */
/*                                                          */
/************************************************************/
int busca_diferente_medio(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom)
{
  int i,j,k;
  int mas_diferente, diferencia, mayor_dif;
  int dif_medio, media_dif;

  mayor_dif = 0;
  mas_diferente = 0;
  for(i=0; i < tam_pob; i++){
    diferencia = 0;
    for(k=0; k<NUM_VAR;k++){
      for(j=0; j < lcrom; j++){
        if(pob[i].vec_genes[k][j] != *(bool *)(representa_genotipo[k] + j))
          diferencia++;
      }
    }
    if (diferencia > mayor_dif){
      mayor_dif = diferencia;
      mas_diferente = i;
    }
  }
  media_dif = (int)mas_diferente/2;
  dif_medio = mas_diferente;
  for(i=0; i < tam_pob; i++){
    diferencia = 0;
    for(k=0; k<NUM_VAR;k++){
      for(j=0; j < lcrom; j++){
        if(pob[i].vec_genes[k][j] != *(bool *)(representa_genotipo[k] + j))
          diferencia++;
      }
    }
    if (diferencia == media_dif){
      dif_medio = i;
      break;
    }
  }
  return dif_medio;
}
/************************************************************/
/*                                                          */
/* funcion: busca_diferente_elite                           */
/* Busca en la poblaci'on la posicion del individuo mas     */
/* diferente a represeta_genotipo  y que este en la elite   */
/*                                                          */
/************************************************************/
int busca_diferente_elite_fija(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom, Tvec_int& orden_pob,
int tam_elite)
{
  int i,j,k,l;
  int mas_diferente, diferencia, mayor_dif;

  // se busca al primer elemento de la elite
/*
  k = 0;
  while((k< tam_pob) && (pob[k].elite == false))
    k++;
*/

  mayor_dif = 0;
  mas_diferente = 0;
/*
  for(i=k; i < tam_pob; i++){
    diferencia = 0;
    for(l=0; l < NUM_VAR; l++){
      for(j=0; j < lcrom; j++){
        if(pob[i].vec_genes[l][j] != *(bool *)(representa_genotipo[l]+j))
          diferencia++;
      }
    }
    if ((diferencia > mayor_dif) && (pob[i].elite == true)){
      mayor_dif = diferencia;
      mas_diferente = i;
    }
  }
*/
  //for(i=k; i < (int)tam_pob/4; i++){
  //for(i=k; i < (int)tam_pob/reduccion_elite; i++){
  //for(i=1; i < (int)tam_pob/reduccion_elite; i++){
  for(i=1; i < tam_elite; i++){
    diferencia = 0;
    for(l=0; l < NUM_VAR; l++){
      for(j=0; j < lcrom; j++){
        if(pob[orden_pob[i]].vec_genes[l][j] != *(bool *)(representa_genotipo[l]+j))
          diferencia++;
      }
    }
    if (diferencia > mayor_dif){
      mayor_dif = diferencia;
      mas_diferente = orden_pob[i];
    }
  }
  return mas_diferente;
}

/************************************************************/
/*                                                          */
/* funcion: busca_diferente_elite                           */
/* Busca en la poblaci'on la posicion del individuo mas     */
/* diferente a represeta_genotipo  y que este en la elite   */
/*                                                          */
/************************************************************/
int busca_diferente_elite(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom, Tvec_int& orden_pob,
int reduccion_elite)
{
  int i,j,k,l;
  int mas_diferente, diferencia, mayor_dif;

  // se busca al primer elemento de la elite
/*
  k = 0;
  while((k< tam_pob) && (pob[k].elite == false))
    k++;
*/

  mayor_dif = 0;
  mas_diferente = 0;
/*
  for(i=k; i < tam_pob; i++){
    diferencia = 0;
    for(l=0; l < NUM_VAR; l++){
      for(j=0; j < lcrom; j++){
        if(pob[i].vec_genes[l][j] != *(bool *)(representa_genotipo[l]+j))
          diferencia++;
      }
    }
    if ((diferencia > mayor_dif) && (pob[i].elite == true)){
      mayor_dif = diferencia;
      mas_diferente = i;
    }
  }
*/
  //for(i=k; i < (int)tam_pob/4; i++){
  //for(i=k; i < (int)tam_pob/reduccion_elite; i++){
  for(i=1; i < (int)tam_pob/reduccion_elite; i++){
    diferencia = 0;
    for(l=0; l < NUM_VAR; l++){
      for(j=0; j < lcrom; j++){
        if(pob[orden_pob[i]].vec_genes[l][j] != *(bool *)(representa_genotipo[l]+j))
          diferencia++;
      }
    }
    if (diferencia > mayor_dif){
      mayor_dif = diferencia;
      mas_diferente = orden_pob[i];
    }
  }
  return mas_diferente;
}

/************************************************************/
/*                                                          */
/* funcion: busca_diferente                               */
/* Busca en la poblaci'on la posicion del individuo mas     */
/* diferente a represeta_genotipo                           */
/*                                                          */
/************************************************************/
void calcula_secuencia_consenso(TPoblacion pob, int tam_pob, char** secuencia_consenso, int lcrom)
{
  int i,j,k;
  int** cuenta_veces;   // cuenta el numeros de 1 en cada gen

  cuenta_veces = new int*[NUM_VAR];
  for(k=0; k<NUM_VAR; k++)
    cuenta_veces[k] = new int[lcrom];

  for(k=0; k<NUM_VAR; k++)
    for(j=0; j < lcrom; j++)
      cuenta_veces[k][j] = 0;
  for(i=0; i < tam_pob; i++){
    for(k=0; k<NUM_VAR; k++){
      for(j=0; j < lcrom; j++){
        if(pob[i].vec_genes[k][j] == 1)
          cuenta_veces[k][j]++;
      }
    }
  }
//cout << "CUENTA VECES: " << endl;
  //for(j=0; j < lcrom; j++){
    //cout << cuenta_veces[j] << " ";
  //}
//cout << endl;
  for(k=0; k<NUM_VAR; k++){
    for(j=0; j < lcrom; j++){
      if(cuenta_veces[k][j] > (int)(tam_pob / 2))
        //*(char *)(secuencia_consenso[k]+j) = (char)1;
        secuencia_consenso[k][j] = (char)1;
      else
        //*(char *)(secuencia_consenso[k]+j) = (char)0;
        secuencia_consenso[k][j] = (char)0;
    }
  }
//cout << "SECUENCIA CONSENSO: " << endl;
  //for(j=0; j < lcrom; j++){
    //if(secuencia_consenso[j] == 1)
      //cout << "1";
    //else if(secuencia_consenso[j] == 0)
      //cout << "0";
  //}
//cout << endl;
  delete [] cuenta_veces;
}

/************************************************************/
/*                                                          */
/* funcion: buscar_cadena_cache                             */
/* Busca una cadena de genes en la cache y si la encuentra  */
/* devuelve su fitness                                      */
/*                                                          */
/************************************************************/
bool buscar_cadena_cache(TCache cache, TIndividuo individuo, double& fitness)
{
  TIter_cache iter_cache;
  long d;  // valor deciaml de la cadena binaria
  //d = bin_dec(genes, LON_CAD);
//cout << "BUSCAR CADENA: d = " << d << endl;
  string cad_genes;
  int i,j;

  //d = (double)bin_dec(genes, LON_CAD);
  for(j=0; j < NUM_VAR; j++){
    for(i=0; i < LON_CAD; i++){
      if (individuo.vec_genes[j][i] == 1)
        cad_genes.append("1");
      else
        cad_genes.append("0");
    }
  }

  //iter_cache = cache.find(d);
  iter_cache = cache.find(cad_genes);
  if (iter_cache == cache.end())
    return false;
  else{
    fitness = (iter_cache->second);
    return true;
  } 
}

/************************************************************/
/*                                                          */
/* funcion: inserta_cadena_cache                            */
/* Inserta una nueva cadena, con su fitness, en la cache    */
/*                                                          */
/************************************************************/
void insertar_cadena_cache(TCache& cache, TIndividuo individuo, double fitness)
{
  int i,j;
  //long d;  // valor deciaml de la cadena binaria
  double d;  // valor deciaml de la cadena binaria
  string cad_genes;

  //d = (double)bin_dec(genes, LON_CAD);
  for(j=0; j < NUM_VAR; j++){
    for(i=0; i < LON_CAD; i++){
      if (individuo.vec_genes[j][i] == 1)
        cad_genes.append("1");
      else
        cad_genes.append("0");
    }
  }

  //cache.insert(TCache::value_type(d, fitness));
  cache.insert(TCache::value_type(cad_genes, fitness));

//PRUEBA
/*
  TIter_cache iter_cache;
  double fitness2;

  //iter_cache = cache.find(d);
  iter_cache = cache.find(cad_genes);

  fitness2 = (iter_cache->second);
  if (fitness != fitness2){
cout << "ERRORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR!" << endl;
cout << "fitness = " << fitness << endl;
cout << "fitness2 = " << fitness2 << endl;
cout << "iter_cache->first = " << iter_cache->first << endl;
  }
*/
}

/************************************************************/
/*                                                          */
/* funcion: calcula entropia                                */
/* Calcula la entropia de la poblacion en la generacion     */
/* en curso                                                 */
/*                                                          */
/************************************************************/
double calcula_entropia(TPoblacion pob, int tam_pob, int lcrom)
{
  int i,j,k;
  double entropia;
  int num_veces;
  float frecuencia;

  entropia = 0;
  for(i=0; i < tam_pob; i++){
    num_veces = 0;
    for(j=0; j < tam_pob; j++){
      //if ((i != j) && (pob[i].aptitud == pob[j].aptitud))
      if (pob[i].aptitud == pob[j].aptitud)
        num_veces++;
    }
    //cout << "num_veces = " << num_veces << endl;
    //entropia = entropia + num_veces * log((double)num_veces) / log((double)2);
    frecuencia = (float)num_veces / (float)tam_pob;
    //entropia = entropia + frecuencia * log((double)frecuencia);
    entropia = entropia + log((double)frecuencia);
  }
  //entropia = -entropia;
  entropia = -entropia / (float)tam_pob;
  //cout << "entropia = " << entropia << endl;
  return entropia;
}

