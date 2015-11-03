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
public class Photometry {
    
    public static double[] UBVRI(double[] lambdaScale, double[] flux) {

        double[][][] filters = FilterSet.filterSet();

        int numLams = lambdaScale.length;
        double[] logFlux = new double[numLams];
        for (int i = 0; i < numLams; i++){
            logFlux[i] = Math.log(flux[i]);
        }
        
        int numCols = 5;  //five band combinations in Johnson-Bessell UxBxBVRI: Ux-Bx, B-V, V-R, V-I, R-I
        double[] colors = new double[numCols];

        int numBands = filters.length;
        int numLambdaFilt = filters[0][0].length;

        double[] bandFlux = new double[numBands];

        //Note: Calibration must be re-done!  May 2015
        // Single-point Johnson UBVRI calibration to Vega:
        // Vega colours computed self-consistntly with GrayFox 1.0 using 
        // Stellar parameters of Castelli, F.; Kurucz, R. L., 1994, A&A, 281, 817
        // Teff = 9550 K, log(g) = 3.95, ([Fe/H] = -0.5 - not directly relevent):
        //Colours w. 14 spectral lines, 28 Jul 2015:
        //double[] vegaColors = {-0.059805, -0.182563, 0.181314, -0.357444, -0.538758};
        //double[] vegaColors = {0.0, 0.0, 0.0, 0.0, 0.0}; //For re-calibrating with raw Vega colours
        // Aug 2015 - with 14-line linelist:
        double[] vegaColors = {0.289244, -0.400324, 0.222397, -0.288568, -0.510965};

        double deltaLam, newY, product;

        for (int ib = 0; ib < numBands; ib++) {

            bandFlux[ib] = 0.0; //initialization

            //wavelength loop is over photometric filter data wavelengths
            for (int il = 1; il < numLambdaFilt; il++) {

                //In this case - interpolate model SED onto wavelength grid of given photometric filter data
                deltaLam = filters[ib][0][il] - filters[ib][0][il - 1];  //nm
                //deltaLam = 1.0e-7 * deltaLam;  //cm

                //hand log flux (row 1) to interpolation routine: 
                
                newY = Interpol.interpol(lambdaScale, logFlux, filters[ib][0][il]);
                // linearize interpolated flux: - fluxes add *linearly*
                newY = Math.exp(newY);
                //System.out.println("Photometry: newFlux: " + newFlux + " filterlamb: " + filters[ib][0][il]);

                product = filters[ib][1][il] * newY;

                //System.out.println("Photometry: filtertrans: " + filters[ib][1][il] + " product: " + product + " deltaLam: " + deltaLam);
                //Rectangular picket integration
                bandFlux[ib] = bandFlux[ib] + (product * deltaLam);
                //System.out.println("Photometry: ib: " + ib + " bandFlux: " + bandFlux[ib]);

            } //il loop - lambdas
            //System.out.println("Photometry: ib: " + ib + " bandFlux: " + bandFlux[ib]);

        }  //ib loop - bands

        double raw;

        // Ux-Bx: 
        raw = 2.5 * Math.log10(bandFlux[1] / bandFlux[0]);
        colors[0] = raw - vegaColors[0];
        System.out.println("U-B: " + colors[0]);

        // B-V:
        raw = 2.5 * Math.log10(bandFlux[3] / bandFlux[2]);
        colors[1] = raw - vegaColors[1];
        System.out.println("B-V: " + colors[1]);

        // V-R:
        raw = 2.5 * Math.log10(bandFlux[4] / bandFlux[3]);
        colors[2] = raw - vegaColors[2];
        System.out.println("V-R: " + colors[2]);

        // V-I:
        raw = 2.5 * Math.log10(bandFlux[5] / bandFlux[3]);
        colors[3] = raw - vegaColors[3];
        System.out.println("V-I: " + colors[3]);

        // R-I:
        raw = 2.5 * Math.log10(bandFlux[5] / bandFlux[4]);
        colors[4] = raw - vegaColors[4];
        System.out.println("R-I: " + colors[4]);

        return colors;

    }
    
    
}
