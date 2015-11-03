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
public class PhysData {

    //public static void physData(double xFrac, double yFrac, double zFrac){
                    //Physical data:
    // EOS quantities
    //Adiabatic gammas
    public static double gammaMono = 5.0 / 3.0; //ideal monatomic gas 
    public static double gammaRad = 4.0 / 3.0; //radiation (photon gas)
    //Constant "K" in polytropic (adiabatic) EOS (P = K*rho^gamma)
    // Value based on conditions at centre of Sun:
    //public static double polyK = 5.134615452370382E13;
    //Mean molecular weight of fully ionized gas of solar composition
    //double muI = 0.62; // amu (C&O 2nd Ed. p. 293)
    //The following assumes Sigma_elements((1+z)/A) ~ 0.5:
    //C&O 2nd Ed., p. 293

    public static double muI(double xFrac, double yFrac, double zFrac) {
        double invMuI = 2.0 * xFrac + 0.75 * yFrac + 0.5 * zFrac;
        return 1.0 / invMuI;
    }

    public static double logGamMono(double gammaMono) {
        return Math.log(gammaMono);
    }

    public static double logGamRad(double gammaRad) {
        return Math.log(gammaRad);
    }

    public static double logMuI(double xFrac, double yFrac, double zFrac) {
        return Math.log(PhysData.muI(xFrac, yFrac, zFrac));
    }

    //public static double logPolyK(){
    //    return Math.log(polyK);
    //}    
    //Kramers oapcity pre-factors (C&O 2nd Ed. p. 249-250):
    // Assume all Gaunt and guillotine factors (t) are 1.0
    // units: cm^2 g^-1 (cgs)
    public static double kap0BF = 4.34e22; // bound-free (photo-ionization)
    public static double kap0FF = 3.68e19; // free-free (bremstrahlung)
    public static double kap0ES = 0.2;  // electron scattering
    public static double kap0Hmin = 7.9e-33 / 0.02;  // H^- b-f

    public static double logKap0BF() {
        return Math.log(kap0BF);
    }

    public static double logKap0FF() {
        return Math.log(kap0FF);
    }

    public static double logKap0ES() {
        return Math.log(kap0ES);
    }

    public static double logKap0Hmin() {
        return Math.log(kap0Hmin);
    }

    //Nuclear E generation data (C&O 2nd Ed. p. 311 - 312):
    // Power law mass power generation rate pre-factors
    //proton-proton (p-p) chain H fusion
    public static double eps0PP = 1.08e-5; // ergs s^-1 cm^3 g^-2
    public static double betaPP = 4.0; //T_6 exponent

    public static double logEps0PP() {
        return Math.log(eps0PP);
    }
    //CNO cycle - H fusion catalyzed by Carbon, Nitrogen and Oxygen 
    public static double eps0CNO = 8.24e-24;  // ergs s^-1 cm^3 g^-2
    public static double betaCNO = 19.9; //T_6 exponent

    public static double logEps0CNO() {
        return Math.log(eps0CNO);
    }

    // threshold for H-fusion
    //C&O 2nd Ed., p. 302
    public static double fusionTemp = 1.0e7; //K 

    //}
}
