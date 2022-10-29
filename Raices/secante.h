#ifndef SECANTE_H
#define SECANTE_H

#include <functional>
#include "APROXIMACION.H"

using std::function;
using raices::raiz;
	
namespace raices
{
	class secante
	{
		public:
			secante(function <double(double)> prmFuncion):atrFuncion(prmFuncion)
			{
			}
			raiz calcularRaiz(double prmXo, double prmX1, double erp, int maxIter)
			{
				int i = 0;
				raiz r;
				while(i <= maxIter)
				{
					double p2 = prmX1 - (atrFuncion(prmX1)*(prmXo-prmX1))/
											(atrFuncion(prmXo)-atrFuncion(prmX1));
					if(i > 0)
					{
						r.agregar(p2);
					}
					double vErp = fabs((p2-prmX1)/p2)*100.0f;
					if(vErp < erp)
					{
						r.setValor(p2);
						return r;
					}
					i = i+1;
					prmXo = prmX1;
					prmX1 = p2;
				}
				return r;
			}
		private:
			function <double(double)> atrFuncion;
		
	};
};

#endif
