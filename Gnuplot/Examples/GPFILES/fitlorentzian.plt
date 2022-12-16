set yrange [*:*]


fitfunction(x) = (fitlambda * fitgamma* fitgamma) / ( ( (x-fitenergy)*(x-fitenergy)) + (fitgamma* fitgamma*(1+fitlambda)*(1+fitlambda)/4) )

fit fitfunction(x) INPUT_FILE using ($1-FermiEnergy):2 via fitlambda, fitgamma, fitenergy


seebeckcoefficient(x) = 2*S0 * (x-fitenergy) * fitfunction(x) / (fitlambda * fitgamma* fitgamma)

ZTK(x) = S0 *S0 * G0 * 4* fitfunction(x) * fitfunction(x)  * fitfunction(x) * (x-fitenergy) * (x-fitenergy)  / ((fitlambda * fitgamma* fitgamma) * (fitlambda * fitgamma* fitgamma))

#print(Scaling_factor_Y)
#print(INPUT_FILE)
#plot INPUT_FILE  u ($1+FermiEnergy):2 with filledcurve y1=0  lt 5 lw 2, INPUT_FILE u ($1+FermiEnergy):2  lt -1 lw 1
#plot fitfunction(x) u 1:2 with points  2
