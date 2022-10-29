#ifndef NEWTONGENERALIZADO_H
#define NEWTONGENERALIZADO_H

#include <functional>

using std::function;

namespace raices
{
	class newtonGeneralizado
	{
	public:
		newtonGeneralizado(function <double(double)> prmFuncion, 
						   function <double(double)> prmDerivada, 
						   function <double(double)> prmSDerivada):
				atrFuncion(prmFuncion), atrDerivada(prmDerivada), 
				atrSDerivada(prmSDerivada)
		{			
		}
		raiz calcularRaiz(double prmPunto, double erp, int maxIter)
		{
			int i = 1;
			raiz r;
			while(i <= maxIter)
			{
				double p = prmPunto - ((atrFuncion(prmPunto)*atrDerivada(prmPunto))/
										(pow(atrDerivada(prmPunto),2.0f)-(atrFuncion(prmPunto)*atrSDerivada(prmPunto))));
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
		function <double(double)> atrSDerivada;
		
	};
};
#endif
