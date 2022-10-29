#ifndef APROXIMACION_H
#define APROXIMACION_H

#include <vector>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;

using std::vector;
namespace raices
{
	struct raiz
	{
		bool encontrada;
		double valor;
		vector<double> aproximaciones;
		int iteraciones;
		vector<double> erpCalculado;
		raiz()
		{
			encontrada = false;
			valor = NAN;
			iteraciones = 0;
		}
	
		void agregar(double prmValor)
		{
			aproximaciones.push_back(prmValor);
			if(iteraciones >= 1)
			{
				erpCalculado.push_back(fabs((aproximaciones[iteraciones-1]-aproximaciones[iteraciones])/aproximaciones[iteraciones])*100.0f);
			}
			else{
				erpCalculado.push_back(NAN);
			}
			iteraciones++;
			
		}
		void setValor(double prmValor)
		{
			valor = prmValor;
			encontrada = true;
		}
		void imprimir()
		{
			if(encontrada)
			{
				for(int i = 0; i < iteraciones; i++)
				{
					cout<<" Iteracion "<<i+1<<" valor: "<<aproximaciones[i]<<" ERP: "<< std::setprecision(10) <<erpCalculado[i] <<endl;
				}
			}
		}
	};
	
};
#endif
