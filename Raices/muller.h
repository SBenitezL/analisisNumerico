#ifndef MULLER_H
#define MULLER_H
#include <functional>
#include "APROXIMACION.H"

using raices::raiz;
using std::isnan;

namespace raices
{
	class muller
	{
		public:
			
			muller(function<double(double)> prmFuncion):atrFuncion(prmFuncion)
			{
			}
			raiz calcularRaiz(double prmX0, double prmX1, double prmX2,
						 double prmERP, int maxIter)
			{
				//Declaración de variables locales
				double D,E,h,p,erp,b;
				//Contador de iteraciones
				int i = 0;
				//Estructura de datos que almacena las aproximaciones
				raiz r;
				//Calculo de las h
				double h1 = prmX1-prmX0, h2 = prmX2-prmX1;
				//Calculo de los deltha
				double d1= (atrFuncion(prmX1)-atrFuncion(prmX0))/h1,
					   d2 = (atrFuncion(prmX2)-atrFuncion(prmX1))/h2;
				//Calculo de a
				double a=(d2-d1)/(h1+h2);
				
				i = 2;
				while (i < maxIter)
				{					
					b = d2+(h2*a);
					D = sqrt((b*b)-(4.0f*atrFuncion(prmX2)*a));
					//Se toma el mayor para tener menos iteraciones.
					if(fabs(b-D) < fabs(b+D))
					{
						E = b+D;
					}else{
						E  = b-D;
					}
					h = (-2.0f*atrFuncion(prmX2))/E;
					p = prmX2 + h;
					r.agregar(p);
					if(isnan(p))
					{
						r.setValor(NAN);
						return r;
					}
					erp = fabs((p-prmX2)/p)*100.0f;
					if(erp <= prmERP)
					{
						r.setValor(p);
						return r;
					}
					prmX0 = prmX1;
					prmX1 = prmX2;
					prmX2 = p;
					h1 = prmX1-prmX0;
					h2 = prmX2-prmX1;
					d1= (atrFuncion(prmX1)-atrFuncion(prmX0))/h1;
					d2 = (atrFuncion(prmX2)-atrFuncion(prmX1))/h2;
					a=(d2-d1)/(h1+h2);
					i++;
				}
			r.setValor(0.01f);
			return r;
			}
		private:
			function<double(double)> atrFuncion;
	};
};


#endif
