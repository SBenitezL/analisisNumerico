#ifndef STEFFENSON_H
#define STEFFENSON_H

#include <functional>
#include "APROXIMACION.H"

namespace raices
{
	class steffenson
	{
		public:	
			steffenson(function<double(double)> prmFuncion):atrFuncion(prmFuncion)
			{
			}
			raiz calcularRaiz(double prmP,double prmERP, int maxIter)
			{
				int i = 1;
				raiz r;
				double p;
				double erp;
				while(i <= maxIter)
				{
					p = prmP -(pow(atrFuncion(prmP),2)/
							(atrFuncion(prmP+atrFuncion(prmP))-atrFuncion(prmP)));
					r.agregar(p);
					erp = fabs((p-prmP)/p)*100;
					if(erp <= prmERP)
					{
						r.setValor(p);
						return r;
					}
					i++;
					prmP = p;
				}
				return r;
			}
		private:
			function<double(double)> atrFuncion;
	};
};
#endif
