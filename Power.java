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
    
    public static double ppChain(double temp, double rho, double xFrac, double zFrac){ //cnoFrac
                
        // H fusion
        
        // Assumes power law approximation for reaction and power rates
        // Assumes screening factor, corrections for PPII and PPIII and
        // higher order corrections are all unity
        // C&O 2nd Ed. p. 311-312
        
        //p-p chain
        
        double logRatePP;
        double logT6 = Math.log(1.0e-6*temp);  //log "T6" temperature
        double logRho = Math.log(rho);
        double logX = Math.log(xFrac);
        
        logRatePP = PhysData.logEps0PP() + logRho + 2.0 * logX 
                + PhysData.betaPP * logT6;
        
        
        return Math.exp(logRatePP);
        
    }
    
    public static double cnoCycle(double temp, double rho, double xFrac, double zFrac){ //cnoFrac
        
        // H fusion
        
        // Assumes power law approximation for reaction and power rates
        // Assumes screening factor, corrections for PPII and PPIII and
        // higher order corrections are all unity
        // C&O 2nd Ed. p. 311-312
        
        // CNO cycle
        // Need value for X_CNO mass fraction!
        double logT6 = Math.log(1.0e-6*temp);  //log "T6" temperature
        double logRho = Math.log(rho);
        double logX = Math.log(xFrac);       
        double logRateCNO;  
        double logXCNO = Math.log(zFrac/2.0);
        
        logRateCNO = PhysData.logEps0CNO() + logRho + logX + logXCNO 
                + PhysData.betaCNO * logT6;
        
        //System.out.println("totEpsilon " + totEpsilon);
        return Math.exp(logRateCNO);
        
    }     
    
}
