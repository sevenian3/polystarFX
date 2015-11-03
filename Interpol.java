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
public class Interpol {
    
    public static double interpol(double[] x, double[] y, double newX) {

        double newY;

        // Bracket newX:
        double x1, x2;
        int p1, p2;
        p1 = 0;
        p2 = 1;
        x1 = x[p1];
        x2 = x[p2];

        for (int i = 1; i < x.length; i++) {
            if (x[i] >= newX) {
                // Found upper bracket
                p2 = i;
                p1 = i - 1;
                x2 = x[p2];
                x1 = x[p1];
                break;
            }
        }

        double step = x2 - x1;

    //Interpolate
        //First order Lagrange formula
        //   newY = y[1][p2] * (newX - x1) / step
        //           + y[1][p1] * (x2 - newX) / step;
        newY = y[p2] * (newX - x1) / step
                + y[p1] * (x2 - newX) / step;

        //System.out.println("Interpol: p1, p2, x1, x2, y1, y2, newX, newY: " + 
        //        p1 + " " + p2 + " " + x1 + " " + x2 + " " + y[1][p1] + " " + y[1][p2] + " " + newX + " " + newY + " ");
        return newY;

    }    
    
}
