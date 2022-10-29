#ifndef BISECCION_H
#define BISECCION_H

#include <functional>
#include "APROXIMACION.H"

using std::function;
using raices::raiz;

namespace raices
{
	class biseccion
	{
		public: 
			biseccion(function <double(double)> prmFuncion):atrFuncion(prmFuncion)
			{
			}
			raiz calcularRaiz(
								double Xa, 
								double Xb, 
								double erp,
								int maxIter)
			{
				raiz r;
				int i=1;
				double Xanterior = (Xa + Xb)/2.0f;
				if(atrFuncion(Xa) * atrFuncion(Xanterior) <0)
				{
					//Si tienen el mismo signo quiere decir que la raíz se 
					//encuentra en el segundo intervalo
					Xb = Xanterior;
				}
				else
				{
					//Esta en el primer intervalo
					Xa = Xanterior;
				}
				while(i<=maxIter)
				{
					double Xnueva = (Xa + Xb)/2.0f;
					double er = fabs((Xnueva - Xanterior)/Xnueva)*100.0f;
					r.agregar(Xnueva);
					if(er <= erp)
					{
						r.setValor(Xnueva);
						return r;
					}
					i = i+1;
					if(atrFuncion(Xa) * atrFuncion(Xnueva) >0)
					{
						//Si tienen el mismo signo quiere decir que la raíz se 
						//encuentra en el segundo intervalo
						Xa = Xnueva;
					}
					else
					{
						//Esta en el primer intervalo
						Xb = Xnueva;
					}
					Xanterior = Xnueva;
				}
				return r;
			}
		private:
			function<double(double)> atrFuncion;
	};
};
#endif
