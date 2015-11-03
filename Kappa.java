/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package polystar;

/**
 *
 * @author Ian
 */
public class Kappa {

    public static double kappaBfFn(double temp, double rho, double xFrac, double zFrac) {

        double kappa;

        double logTemp = Math.log(temp);
        double logRho = Math.log(rho);

        //Kramers-type oapcity laws
        //Assumes all Gaunt and guillotine (cut-off) factors are unity
        //C&O 2nd Ed. p. 249-250
        //Combined B-F Gaunt and Guillotine (cut-off) cutoff factors:
        double logGtBF = Math.log(0.708) + 0.2 * (Math.log(rho) + Math.log(1.0 + xFrac));

        //bounf-free (b-f) contribution (photo-ionization)
        double logKapBF = PhysData.logKap0BF() - logGtBF + logRho - 3.5 * logTemp
                + Math.log(zFrac) + Math.log(1.0 + xFrac);
        double kapBF = Math.exp(logKapBF);

        return kapBF;

    }

    public static double kappaFfFn(double temp, double rho, double xFrac, double zFrac) {

        double logTemp = Math.log(temp);
        double logRho = Math.log(rho);

        //Kramers-type oapcity laws
        //Assumes all Gaunt and guillotine (cut-off) factors are unity
        //C&O 2nd Ed. p. 249-250
        //free-free (f-f) contribution (bremstrahlung, braking radiation)
        double logKapFF = PhysData.logKap0FF() + logRho - 3.5 * logTemp
                + Math.log(1.0 - zFrac) + Math.log(1.0 + xFrac);
        double kapFF = Math.exp(logKapFF);

        //System.out.println("kappa " + kappa);
        return kapFF;

    }

    public static double kappaEsFn(double temp, double rho, double xFrac, double zFrac) {

        double logTemp = Math.log(temp);
        double logRho = Math.log(rho);

        //Kramers-type oapcity laws
        //Assumes all Gaunt and guillotine (cut-off) factors are unity
        //C&O 2nd Ed. p. 249-250        
        //Thomson scattering from free electrons:
        double logKapES = PhysData.logKap0ES() + Math.log(1.0 + xFrac);
        double kapES = Math.exp(logKapES);

        return kapES;

    }

    public static double kappaHminFn(double temp, double rho, double xFrac, double zFrac) {

        double logTemp = Math.log(temp);
        double logRho = Math.log(rho);
        //initializations:
        double logKapHmin = -99.0;
        double kapHmin = 0.0;

        //Kramers-type oapcity laws
        //Assumes all Gaunt and guillotine (cut-off) factors are unity
        //C&O 2nd Ed. p. 249-250   
        // Hminus b-f:
        // Hminus opacity is DANGEROUS: T^9 !! 
        if ((temp > 3000.0) && (temp < 6000.0)
                && (rho > 1.0e-13) && (rho < 1.0e-8)
                && (zFrac > 0.001) && (zFrac < 0.03)) {
            logKapHmin = PhysData.logKap0Hmin() + 0.5 * logRho + 9.0 * logTemp
                    + Math.log(zFrac);
            kapHmin = Math.exp(logKapHmin);
        } else {
            kapHmin = 0.0; //initialization
        }
        
        return kapHmin;
    }

}
