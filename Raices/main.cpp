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
#include "STEFFENSON.H"
#include "muller.h"

#include <vector>

using std::cout;
using std::endl;
using std::setprecision;

using raices::biseccion;
using raices::reglaFalsa;
using raices::newtonRaphson;
using raices::secante;
using raices::newtonGeneralizado;
using raices::steffenson;
using raices::muller;

using std::string;

using std::vector;
using raices::raiz;

struct strParametros
{
	vector<double> valoresIniciales;
	double erp;
	function <double(double)> funcion;
	function <double(double)> derivadaf;
	function <double(double)> sDerivadaf;
	int maxIter;
	strParametros(vector<double> prmValores,function<double(double)> prmFuncion, 
				  function<double(double)> prmDerivada, function <double(double)> prmSDerivada, double prmErp,
				  int prmIter): valoresIniciales(prmValores),funcion(prmFuncion),
		derivadaf(prmDerivada), sDerivadaf(prmSDerivada),maxIter(prmIter)
		
	{
		erp = prmErp;
	}
};
vector<raiz> pruebas(strParametros prmParametros)
{
	biseccion bs(prmParametros.funcion);
	reglaFalsa rF(prmParametros.funcion);
	newtonRaphson nR(prmParametros.funcion,prmParametros.derivadaf);
	secante sct(prmParametros.funcion);
	newtonGeneralizado ng(prmParametros.funcion, prmParametros.derivadaf,prmParametros.sDerivadaf);
	steffenson stf(prmParametros.funcion);
	
	vector<raiz> resultados;
	resultados.push_back(bs.calcularRaiz(prmParametros.valoresIniciales[0],
								   prmParametros.valoresIniciales[1],
								   prmParametros.erp,prmParametros.maxIter));
	
	resultados.push_back(rF.calcularRaiz(prmParametros.valoresIniciales[0],
								   prmParametros.valoresIniciales[1],
								   prmParametros.erp,prmParametros.maxIter));
	resultados.push_back(nR.calcularRaiz(prmParametros.valoresIniciales[2],
								   prmParametros.erp,prmParametros.maxIter));
	resultados.push_back(sct.calcularRaiz(prmParametros.valoresIniciales[3],
								   prmParametros.valoresIniciales[4],
								   prmParametros.erp,prmParametros.maxIter));
	resultados.push_back(ng.calcularRaiz(prmParametros.valoresIniciales[2],
										prmParametros.erp,prmParametros.maxIter));
	resultados.push_back(stf.calcularRaiz(prmParametros.valoresIniciales[2],
										  prmParametros.erp,prmParametros.maxIter));
	return resultados;
}
int main (int argc, char *argv[]) {
	string metodos[] = {"Biseccion","Regla Falsa","Newton Raphson","Secante","Newton Generalizado","Steffenson"};
	strParametros varParametros({-0.5f,1.0f,1.5f,0.0f,1.5f},
		[](double x)-> double
		{
			return pow(x,3.0f)+(4.0f*pow(x,2.0f))-10.0f;
		},
		[](double x)-> double
		{
			return (3.0f*pow(x,2.0f))+(8.0f*x);
		},
		[](double x)-> double
		{
			return (6.0f*x)+8.0f;
		},0.5f,1000);
	vector<raiz> resultados = pruebas(varParametros);
	for( int i = 0; i < 6; i++)
	{
		cout <<metodos[i]<<": "<< endl<<endl;
		resultados[i].imprimir();
		cout <<"-----------------------------------------------------"<<endl;
	}
	muller ml (
			   [](double x)-> double
			   {
				   return pow(x,3.0f)+(4.0f*pow(x,2.0f))-10.0f;
			   }
			   );
	cout<<"Muller: "<<endl;
	raiz raizML = ml.calcularRaiz(1.5f,2.0f,2.5f, 1.0f,100);
	raizML.imprimir();
}

