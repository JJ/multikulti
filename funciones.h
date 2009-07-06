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

#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <string>    
#include <vector.h>
#include <map.h>
#include "pvm3.h"

using namespace std;

int NUM_VAR = 15; // para problema MMDP(15 var, 90 bits)
const int LON_CAD = 6;  // numero de bits de la cadena para MMDP
//int NUM_VAR = 10; // para problema MMDP(10 var, 60 bits)
//const int LON_CAD = 6;  // numero de bits de la cadena para MMDP
//int NUM_VAR = 1; // para problema P-peaks y wP-Peaks
//const int LON_CAD = 100;  // numero de bits de la cadena para PPeaks
//const int NUM_PICOS = 100;  // numero de bits de la cadena para PPeaks
//int NUM_VAR = 1; // para problema P-peaks y wP-Peaks
//const int LON_CAD = 30;  // numero de bits de la cadena para wPPeaks
//const int NUM_PICOS = 10;  // numero de bits de la cadena para wPPeaks
//int NUM_VAR = 1; // para problema P-peaks y wP-Peaks
//const int LON_CAD = 128;  // numero de bits de la cadena para PPeaks(2)
//const int NUM_PICOS = 128;  // numero de bits de la cadena para PPeaks(2)

typedef vector < bool > TGenes;   // cadena de bits (genotipo)
typedef vector< TGenes > Tvec_genes; // vector de genotipos de cada variable
typedef vector< float > Tvec_var;    // vector de fenotipos de cada variable
typedef vector< int > Tvec_int;    // vector para la longitudes de los genes
typedef struct{
          Tvec_genes vec_genes;
          Tvec_var vec_x;          // fenotipo
          double aptitud;    // funcion de evaluacion
          double puntuacion; // puntuacion relativa : fitness / sumfitness
          double punt_acu;   // puntuacion acumulada para sorteos
          bool elite;       // indica si pertenece a la elite
        }TIndividuo;


typedef vector <TIndividuo> TPoblacion;

/*********************************************************************/
/* clase para representar una cadena de bits que implemente la cache */
/*********************************************************************/
//typedef map< long, double, less<long> > TCache;  // cache para almacenar las fitness ya calculadas 
//typedef map< long, double, less<long> >::const_iterator TIter_cache;  // iterador
typedef map< string, double > TCache;  // cache para almacenar las fitness ya calculadas 
typedef map< string, double >::const_iterator TIter_cache;  // iterador

// funciones de uso de la cache
bool buscar_cadena_cache(TCache cache, TIndividuo individuo, double& fitness);
void insertar_cadena_cache(TCache& cache, TIndividuo individuo, double fitness);

/************************/
/* Funciones del modulo */
/************************/
TPoblacion poblacion_inicial(int tam_pob, int lcrom, float x_min, float x_max);
TIndividuo genera_indiv(int lcrom, float x_min, float x_max);
void seleccion(TPoblacion& pob, int tam_pob);
double adaptacion(TIndividuo& individuo, int lcrom, float x_min, float x_max);
double decod(TGenes genes, int lcrom, float x_min, float x_max); 
int bin_dec(TGenes genes, int lcrom); 

void reproduccion(TPoblacion& pob, int tam_pob, int lcrom, float prob_cruce, 
                  float x_min, float x_max);
void cruce(TIndividuo padre1, TIndividuo padre2, TIndividuo& hijo1, TIndividuo& hijo2,
           int lcrom, int punto_cruce, float x_min, float x_max);  
void cruce_bipunto(TIndividuo padre1, TIndividuo padre2, TIndividuo& hijo1, TIndividuo& hijo2,
           int lcrom, int punto_cruce1, int punto_cruce2, float x_min, float x_max);
void mutacion(TPoblacion& pob, int tam_pob, int lcrom, float prob_mut, float x_min, float x_max);
void evaluacion(TPoblacion& pob, int tam_pob, int& pos_mejor, float& sumaptitud, Tvec_int& orden_pob);
int seleccion_emigrante_mejor(TPoblacion pob, int tam_pob, float sumaptitud);
int alea( int min_valor, int max_valor );
int busca_diferente(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom);
int busca_diferente_medio(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom);
void calcula_secuencia_consenso(TPoblacion pob, int tam_pob, char** secuencia_consenso, int lcrom);
//int busca_diferente_elite(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom);
int busca_diferente_elite(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom, 
            Tvec_int& orden_pob, int reduccion_elite);
int busca_diferente_elite_fija(TPoblacion pob, int tam_pob, char** representa_genotipo, int lcrom, 
            Tvec_int& orden_pob, int tam_elite);
double calcula_entropia(TPoblacion pob, int tam_pob, int lcrom);


