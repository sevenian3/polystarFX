/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package polystar;

/**
 *
 * @author Ian
 *
 * Use Newton-Raphson root finding to find the temperature for which the total
 * of the gas and radiation pressure equals the given pressure:
 *
 */
public class NRtemp {

    public static double getTemp(double press, double rho, double xFrac, double yFrac, double zFrac, double lastTemp) {

        double tolerance = 1.0e-3;  // relative error tolerance for root
        int maxIter = 1000; //maximum number of N-R iterations

        double rootTemp = 1.0;  //initialization

        double cGas = Useful.k * rho / PhysData.muI(xFrac, yFrac, zFrac) / Useful.amu;
        double cRad = Useful.aStef() / 3.0;

        //establish initial guess ("x_0") be averaging the temperature we would have 
        // if the pressure were due entirely to gas and entirely to radiation:
        double temp0 = lastTemp; //initialization
        if (lastTemp == 0.0) {
            double tGas0 = press / cGas;
            double tPress04 = press / cRad;
            double tPress0 = Math.pow(tPress04, 0.25);
            temp0 = 0.5 * (tGas0 + tPress0);
        } else {
            temp0 = lastTemp;
        }

        double y0 = (cRad * Math.pow(temp0, 4.0)) + (cGas * temp0) - press;
        double yPrime0, y1;
        double temp1 = temp0; //initialization of improved guess

        int i;
        for (i = 0; i < maxIter; i++) {

            if (Math.abs(y0) <= tolerance) {
                break;
            } else {

                yPrime0 = (4.0 * cRad * Math.pow(temp0, 3.0)) + cGas;
                temp1 = temp0 - (y0 / yPrime0);

                y1 = (cRad * Math.pow(temp1, 4.0)) + (cGas * temp1) - press;

                temp0 = temp1;
                y0 = y1;

            }
        }

        rootTemp = temp1;

        i--;
        //System.out.println("Number of N-R iterations: " + i);

        return rootTemp;
    }

}
