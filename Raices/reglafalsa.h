#ifndef REGLAFALSA_H
#define REGLAFALSA_H

#include <functional>

using std::function;

namespace raices
{
	class reglaFalsa
	{
		public:
			reglaFalsa (function <double(double)> prmFuncion):atrFuncion(prmFuncion)
			{
			}
			raiz calcularRaiz(double Xi, double Xs, double erp, int maxIter)
			{
				raiz r;
				double Xanterior = Xs - (atrFuncion(Xs)*(Xi-Xs))/(atrFuncion(Xi)-atrFuncion(Xs));
				if(atrFuncion(Xi) * atrFuncion(Xanterior) < 0)
				{
					Xs = Xanterior;
				}
				else
				{
					Xi = Xanterior;
				}
				int i = 1;
				while(i < maxIter)
				{
					double Xnueva = Xs - (atrFuncion(Xs)*(Xi-Xs))/(atrFuncion(Xi)-atrFuncion(Xs));
					double er = fabs((Xnueva-Xanterior)/Xnueva) *100.0f;
					r.agregar(Xnueva);
					if(er <= erp)
					{
						r.setValor(Xnueva);
						return r;
					}
					if(atrFuncion(Xi) * atrFuncion(Xnueva) < 0)
					{
						Xs = Xnueva;
					}
					else
					{
						Xi = Xnueva;
					}
					i = i+1;
					Xanterior = Xnueva;
				}
				return r;
			}
		private:
			function <double(double)> atrFuncion;
	};
};
#endif
