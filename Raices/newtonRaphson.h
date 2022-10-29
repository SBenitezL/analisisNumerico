#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H
 
#include <functional>
#include<iostream>
#include "aproximacion.h"


using std::cout;
using std::endl;
using std::function;
using raices::raiz;

namespace raices
{
	class newtonRaphson
	{
		public:
			newtonRaphson(
						  function<double(double)> prmFuncion,
						  function<double(double)> prmDFuncion
						  ): atrFuncion(prmFuncion), atrDerivada(prmDFuncion)
			{
			}
			raiz calcularRaiz(double prmPunto, double erp, int maxIter)
			{
				int i = 1;
				raiz r;
				while(i <= maxIter)
				{
					double p = prmPunto - (atrFuncion(prmPunto)/atrDerivada(prmPunto));
					r.agregar(p);
					double vErp = fabs((p-prmPunto)/p)*100.0f;
					if(vErp < erp)
					{
						r.setValor(p);
						return r;
					}
					i = i+1;
					prmPunto = p;
				}
				return r;
			}
		private:
		
			function <double(double)> atrFuncion;
			function <double(double)> atrDerivada;
			
	};
	
};
#endif
