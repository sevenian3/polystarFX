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
public class Power {

    public static double ppChain(double temp, double rho, double xFrac, double zFrac) { //cnoFrac

        // H fusion
        // Assumes power law approximation for reaction and power rates
        // Assumes screening factor, corrections for PPII and PPIII and
        // higher order corrections are all unity
        // C&O 2nd Ed. p. 311-312
        //p-p chain
        double logRatePP;
        double T6 = 1.0e-6 * temp;
        double logT6 = Math.log(T6);  //log "T6" temperature
        double logRho = Math.log(rho);
        double logX = Math.log(xFrac);

        //logRatePP = PhysData.logEps0PP() + logRho + 2.0 * logX
         //       + PhysData.betaPP * logT6;
        //
        //More realistic non-power law treatment:

        double fpp = 1.0; //For now...   p-p chain screening factor
        //double psipp = 1.0;  //For now...  correction factor for simultaneous occurence of PPI, PPII & PPIII chains
        //double Cpp = 1.0; //For now... higher order correction terms
        //Notes: from C&) statstar.Physics.Nuclear
        //       PP chains (see Hansen and Kawaler, Eq. 6.65, 6.73, and 6.74)
        //       psipp = 1 + 1.412E8*(1/X - 1)*EXP(-49.98*T6**(-onethird))
        double xTerm = (1.0 / xFrac - 1.0);
        double logPsipp = Math.log(1.412e8) + Math.log(xTerm) - 49.98 * Math.pow(T6, -0.333333);
        double psipp = 1.0 + Math.exp(logPsipp);
        //Cpp = 1 + 0.0123*T6**onethird + 0.0109*T6**twothirds + 0.000938*T6
        double logTerm1 = Math.log(0.0123) + (1.0 / 3.0) * Math.log(T6);
        double logTerm2 = Math.log(0.0109) + (2.0 / 3.0) * Math.log(T6);
        double logTerm3 = Math.log(0.000938) + Math.log(T6);
        double Cpp = 1.0 + Math.exp(logTerm1) + Math.exp(logTerm2) + Math.exp(logTerm3);
        //
        logRatePP = PhysData.logEps2PP() + logRho + 2.0 * logX
                - (2.0 / 3.0) * logT6 - 33.80 * Math.pow(T6, -0.333333)
                + Math.log(fpp) + Math.log(psipp) + Math.log(Cpp);

        return Math.exp(logRatePP);

    }

    public static double cnoCycle(double temp, double rho, double xFrac, double zFrac) { 

        // H fusion
        // Assumes power law approximation for reaction and power rates
        // Assumes screening factor, corrections for PPII and PPIII and
        // higher order corrections are all unity
        // C&O 2nd Ed. p. 311-312
        // CNO cycle
        // Need value for X_CNO mass fraction!
        double T6 = 1.0e-6 * temp;
        double logT6 = Math.log(T6);  //log "T6" temperature
        double logRho = Math.log(rho);
        double logX = Math.log(xFrac);
        double logRateCNO;
        double logXCNO = Math.log(zFrac / 2.0);

        //logRateCNO = PhysData.logEps0CNO() + logRho + logX + logXCNO
        //       + PhysData.betaCNO * logT6;
//
        //More realistic non-power law treatment:
        //double Ccno = 1.0; //For now... higher order correction terms
        //Notes: from C&) statstar.Physics.Nuclear
        //CNO cycle (Kippenhahn and Weigert, Eq. 18.65)
        //CCNO = 1 + 0.0027*T6**onethird - 0.00778*T6**twothirds - 0.000149*T6         
        double logTerm1 = Math.log(0.0027) + (1.0 / 3.0) * Math.log(T6);
        double logTerm2 = Math.log(0.00778) + (2.0 / 3.0) * Math.log(T6);
        double logTerm3 = Math.log(0.000149) + Math.log(T6);
        double Ccno = 1.0 + Math.exp(logTerm1) - Math.exp(logTerm2) - Math.exp(logTerm3);
        //
        logRateCNO = PhysData.logEps2CNO() + logRho + logX + logXCNO
                - (2.0 / 3.0) * logT6 - 152.28 * Math.pow(T6, -0.333333)
                + Math.log(Ccno);

        //System.out.println("totEpsilon " + totEpsilon);
        return Math.exp(logRateCNO);

    }

    public static double TaProcess(double temp, double rho, double yFrac) { 

        // He fusion
        // Assumes screening factor, corrections for PPII and PPIII and
        // higher order corrections are all unity
        // C&O 2nd Ed. p. 311-312
        // Triple alpha process
        // Need value for X_CNO mass fraction!
        double T8 = 1.0e-8 * temp;
        double logT8 = Math.log(T8);  //log "T6" temperature
        double logRho = Math.log(rho);
        double logY = Math.log(yFrac);
        double logRateTa;

        //More realistic non-power law treatment:
        double fTa = 1.0; //For now... Triple alpha process screening factor
        logRateTa = PhysData.logEps2Ta() + 2.0 * logRho + 3.0 * logY
                - 3.0 * logT8 - (44.027 / T8);
           // + Math.log(fpp) + Math.log(psipp) + Math.log(Cpp);       

        //System.out.println("totEpsilon " + totEpsilon);
        return Math.exp(logRateTa);

        //Notes: from C&) statstar.Physics.Nuclear
        //CNO cycle (Kippenhahn and Weigert, Eq. 18.65)
        //CCNO = 1 + 0.0027*T6**onethird - 0.00778*T6**twothirds - 0.000149*T6         
    }

}
