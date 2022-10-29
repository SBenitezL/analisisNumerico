#include<iostream>
#include<cmath>
#include<iomanip>
#include<string>
#include "biseccion.h"
#include "reglafalsa.h"
#include "NEWTONRAPHSON.H"
#include "aproximacion.h"
#include "secante.h"
#include "NEWTONGENERALIZADO.H"
#include "STEFFENSEN.H"
#include "muller.h"

#include <vector>

using std::cout;
using std::endl;
using std::setprecision;
using std::isnan;

using raices::biseccion;
using raices::reglaFalsa;
using raices::newtonRaphson;
using raices::secante;
using raices::newtonGeneralizado;
using raices::steffensen;
using raices::muller;

using std::string;

using std::vector;
using raices::raiz;

struct strParametros
	
#define vERP pow(10.0f,-9.0f)
{
	vector<double> valoresIniciales;
	/*valores iniciales
	*0 Regla falsa
	*1 Regla falsa
	*2 Secante
	*3 Secante
	*4 Muller
	*5 Muller
	*6 Muller
	*7 Steffenson
	*/
	vector<function<double(double)>> funciones;
	/* Funciones
	* 0 Funcion
	* 1 Primera Derivada
	* 2 Segunda Derivada
	* ...
	* n Enésima Derivada
	*/
	double erp;
	int maxIter;
	strParametros(vector<double> prmValores,vector<function<double(double)>> prmFunciones, 
				double prmErp,int prmIter): 
		valoresIniciales(prmValores),funciones(prmFunciones),maxIter(prmIter)
		
	{
		erp = prmErp;
	}
};
vector<raiz> pruebas(strParametros prmParametros)
{
	reglaFalsa rF(prmParametros.funciones[0]);
	secante sct(prmParametros.funciones[0]);
	muller ml(prmParametros.funciones[0]);
	steffensen stf(prmParametros.funciones[0]);
	
	vector<raiz> resultados;
	
	resultados.push_back(rF.calcularRaiz(prmParametros.valoresIniciales[1],prmParametros.valoresIniciales[0],prmParametros.erp,prmParametros.maxIter));
	
	resultados.push_back(sct.calcularRaiz(prmParametros.valoresIniciales[2],prmParametros.valoresIniciales[3],prmParametros.erp,prmParametros.maxIter));
	resultados.push_back(ml.calcularRaiz(prmParametros.valoresIniciales[6],prmParametros.valoresIniciales[5],prmParametros.valoresIniciales[4],prmParametros.erp,prmParametros.maxIter));
	resultados.push_back(stf.calcularRaiz(prmParametros.valoresIniciales[7],prmParametros.erp,prmParametros.maxIter));
	
	return resultados;
}
int main (int argc, char *argv[]) {
	string metodos[] = {"Regla Falsa","Secante","Muller","Steffensen"};
	strParametros varParametros({-16.0f,-14.0f,-20.0f,-10.0f,-13.5f,-14.0f,-14.5f,-14.0f},
		{[](double x)-> double{return 12.2f*(1-exp(-0.04f*x))+5.5f*exp(-0.04f*x);}},
		vERP,1000); 
	vector<raiz> resultados = pruebas(varParametros);
	cout<<endl;
	for( int i = 0; i < 4; i++)
	{
		cout<<" " <<metodos[i]<<": "<< endl<<endl;
		resultados[i].setValor(-15.0f);
		resultados[i].imprimir();
		cout <<"---------------------------------------------------------"<<endl;
	}
}

