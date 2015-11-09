/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package polystar;

import java.text.DecimalFormat;
import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.PasswordField;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.Stage;

/**
 *
 * @author Ian
 */
public class PolyStar extends Application {

    @Override
    public void start(Stage primaryStage) {

        primaryStage.setTitle("PolyStar 1.0 ");

        GridPane grid = new GridPane();
        //grid.setAlignment(Pos.CENTER);
        grid.setHgap(6);
        grid.setVgap(6);
        grid.setPadding(new Insets(15, 15, 15, 15));
        //grid.setAlignment(Pos.TOP_RIGHT);

        Text scenetitle = new Text("PolyStar 1.0 (Hover for tool tips)");
        scenetitle.setFont(Font.font("Tahoma", FontWeight.NORMAL, 16));
        grid.add(scenetitle, 0, 0, 2, 1);

        Label mainLbl = new Label("PolyStar polytropic stellar interior structure modeling");
        mainLbl.setFont(Font.font("Tahoma", FontWeight.NORMAL, 11));
        grid.add(mainLbl, 2, 0);

        ImageView imageLogo;
        imageLogo = new ImageView(
                new Image(PolyStar.class.getResourceAsStream("graphics/SMULOGO2.png")));
        grid.add(imageLogo, 4, 0, 2, 1);

        // Model atmosphere parameters:
        Label atmosLbl = new Label("Model star parameters:");
        atmosLbl.setFont(Font.font("Tahoma", FontWeight.NORMAL, 12));
        grid.add(atmosLbl, 0, 1);
        Tooltip atmosTip = new Tooltip();
        atmosTip.setText(
                "Default parameters are for the Sun's central density, \r\n"
                + " rho_c = 162 g/cm^3, n=1.5"
        );
        atmosLbl.setTooltip(atmosTip);

        //central density input - solar units
        Label rhoCLbl = new Label("rho_c (0.1 - 10.0 x rho_C_Sun)");
        grid.add(rhoCLbl, 0, 2);
        Tooltip rhoCTip = new Tooltip();
        rhoCTip.setText(
                "Central mass density relative to that of Sun \r\n"
                + " rho_C_Sun = 162 g cm^-3"
        );
        rhoCLbl.setTooltip(rhoCTip);

        TextField rhoCSolIn = new TextField("1.0");
        grid.add(rhoCSolIn, 1, 2);

        Label trialLbl = new Label(" ");
        trialLbl.setFont(Font.font("Tahoma", FontWeight.NORMAL, 12));
        grid.add(trialLbl, 2, 1);

        //Polytropic index, n - dimensionless:
        Label indexLbl = new Label("n (1.0 - 3.5)");
        grid.add(indexLbl, 0, 3);
        Tooltip indexTip = new Tooltip();
        indexTip.setText(
                "Dimensionless polytropic index.  \r\n"
                + " = Adiabatic gamma = (n+1)/n \r\n"
                + " n=1.5 for ideal monatomic gas (best for Sun) \r\n"
                + " n=3.0 for photon gas - Eddington standard model"
        );
        indexLbl.setTooltip(indexTip);

        TextField indexIn = new TextField("1.5");
        grid.add(indexIn, 1, 3);
        /*
         // Radiation pressure fraction, beta:
         Label betaLbl = new Label("Beta (0 - 1)");
         grid.add(betaLbl, 0, 4);
         Tooltip betaTip = new Tooltip();
         betaTip.setText(
         "Fractional contribution of gas pressure, P_Gas\r\n"
         + "to total pressure (ie. P_Gas = Beta * P = P - P_Rad)"
         );
         betaLbl.setTooltip(betaTip);

         TextField betaIn = new TextField("0.0");
         grid.add(betaIn, 1, 4);
         */

        //He mass fraction (Y)
        Label yFracLbl = new Label("He mass fraction, Y (0.2 - 0.5)");
        grid.add(yFracLbl, 2, 2);
        Tooltip yFracTip = new Tooltip();
        yFracTip.setText(
                "Helium mass fraction, Y (0.2 - 0.5)"
        );
        yFracLbl.setTooltip(yFracTip);

        TextField yFracIn = new TextField("0.28");
        grid.add(yFracIn, 3, 2);

        //metal mass fraction (Z)
        Label zFracLbl = new Label("Metal mass fraction, Z (0.0001 - 0.05)");
        grid.add(zFracLbl, 2, 3);
        Tooltip zFracTip = new Tooltip();
        zFracTip.setText(
                "Mass fraction of elements of z > 2, Z (0.0001 - 0.05)"
        );
        zFracLbl.setTooltip(zFracTip);

        TextField zFracIn = new TextField("0.02");
        grid.add(zFracIn, 3, 3);

        //
        //
        Button btn = new Button("Model");
        HBox hbBtn = new HBox(10);
        hbBtn.setAlignment(Pos.BOTTOM_RIGHT);
        hbBtn.getChildren().add(btn);
        grid.add(hbBtn, 3, 5);

        final Text actiontarget = new Text();
        grid.add(actiontarget, 0, 8);

        btn.setOnAction(new EventHandler<ActionEvent>() {

            @Override
            public void handle(ActionEvent event) {
                actiontarget.setFill(Color.FIREBRICK);
                //actiontarget.setText("Output here");
                actiontarget.setFont(Font.font("Tahoma", FontWeight.NORMAL, 14));

                String rhoCSolStr = rhoCSolIn.getText();
                String indexStr = indexIn.getText();
                //String betaStr = betaIn.getText();
                String yFracStr = yFracIn.getText();
                String zFracStr = zFracIn.getText();

                // If *any* model atmosphere input string is empty, set rhoCSolStr to empty to simplify 
                // conditional logic below.
                if (indexStr == null || indexStr.isEmpty()) {
                    rhoCSolStr = null;
                }
                //if (betaStr == null || betaStr.isEmpty()) {
                //    rhoCSolStr = null;
                //}
                if (yFracStr == null || yFracStr.isEmpty()) {
                    rhoCSolStr = null;
                }
                if (zFracStr == null || zFracStr.isEmpty()) {
                    rhoCSolStr = null;
                }

                if ((rhoCSolStr != null && !rhoCSolStr.isEmpty())) {

                    //Argument 1: Stellar mass, M, in solar masses
                    double rhoCSol = (Double.valueOf(rhoCSolStr)).doubleValue();

                    //Argument 2: He mass fraftion, Y
                    double yFrac = (Double.valueOf(yFracStr)).doubleValue();

                    //Argument 3: Metal mass fraftion, Z
                    double zFrac = (Double.valueOf(zFracStr)).doubleValue();

                    // Argument 4: Effective temperature, Teff, in K:
                    double index = (Double.valueOf(indexStr)).doubleValue();

                    //// Argument 5: Logarithmic surface gravity, g, in cm/s/s:
                    // double beta = (Double.valueOf(betaStr)).doubleValue();
                    // Sanity checks:
                    if (rhoCSol < 0.1) {
                        rhoCSol = 0.1;
                        rhoCSolStr = "0.1";
                    }
                    if (rhoCSol > 10.0) {
                        rhoCSol = 10.0;
                        rhoCSolStr = "10.0";
                    }

                    if (yFrac < 0.1) {
                        yFrac = 0.1;
                        yFracStr = "0.1";
                    }
                    if (yFrac > 0.5) {
                        yFrac = 0.5;
                        yFracStr = "0.5";
                    }

                    if (zFrac < 0.0002) {
                        zFrac = 0.0002;
                        zFracStr = "0.0002";
                    }
                    if (zFrac > 0.05) {
                        zFrac = 0.05;
                        zFracStr = "0.05";
                    }

                    if (index < 1.0) {
                        index = 1.0;
                        indexStr = "1.0";
                    }
                    if (index > 3.5) {
                        index = 3.5;
                        indexStr = "3.5";
                    }

                    //if (beta < 0.0) {
                    //    beta = 0.0;
                    //    betaStr = "0.0";
                    //}
                    //if (beta > 1.0) {
                    //    beta = 1.0;
                    //    betaStr = "1.0";
                    // }
                    // All code after this line
// Solar parameters:
                    double teffSun = 5778.0;
                    double log10gSun = 4.44;
                    double gravSun = Math.pow(10.0, log10gSun);
//Solar units:
                    double massSun = 1.0;
                    double radiusSun = 1.0;
                    double rhoCSun = 162.2;  // g cm^-3
                    //double massStar = 1.0; //solar masses // test

                    //Composition by mass fraction - needed for opacity approximations
                    //   and interior structure
                    double massXSun = 0.70; //Hydrogen
                    double massYSun = 0.28; //Helium
                    double massZSun = 0.02; // "metals"   
                    // log_10 num density H in surface layer:
                    double log10NH = 17.0;

                    double log10E = Math.log10(Math.E); // convert log_e to log_10
                    double logE10 = Math.log(10.0);     // convert log_10 to log_e

                    double xFrac = 1.0 - yFrac - zFrac;
                    double logXFrac = Math.log(xFrac);
                    double logYFrac = Math.log(yFrac);

                    double rhoC = rhoCSol * rhoCSun; //Sun's central mass density in cgs units

                    //Set up all the special Lane-Emden equation variables:
                    //Adiabatic gamma consistent with input polytropic index
                    double gammaPoly = (index + 1.0) / index;
                    //For convection criterion: 
                    double gamThing = PhysData.gammaMono / (PhysData.gammaMono - 1.0);

                    //Initial Kay and lambda parameters from Lame-Embden equations:
                    /* Doesn't work
                     double Kay, KayTerm1, KayTerm2;
                     KayTerm1 = 3.0 * (1.0 - beta) / Useful.aStef();
                     KayTerm2 = Useful.k / beta / PhysData.muI(xFrac, yFrac, zFrac) / Useful.amu;
                     Kay = Math.pow(KayTerm1, 0.3333) * Math.pow(KayTerm2, 4.0 / 3.0);
                     */
                    double Kay = 5.0e13;  //hard wire to value at centre of Sun

                    double lambda, lamTerm, indxExp;
                    indxExp = (1.0 - index) / index;
                    lamTerm = (index + 1.0) * Kay * Math.pow(rhoC, indxExp) / 4.0 / Math.PI / Useful.GConst;
                    lambda = Math.pow(lamTerm, 0.5);

                    System.out.println("Kay " + Kay + " lambda " + lambda);

                    //Dimensionless Lane-Emden Equation variables, xi and D_n(xi)
                    //xi is the independenet variable
                    //y(xi) is a helper function to separate the 2nd order Lane-Emden equation
                    //into two coupled 1st order equations
                    int maxNumDeps = 1000;
                    double[] xi = new double[maxNumDeps];
                    double[] DFunc = new double[maxNumDeps];
                    double[] yFunc = new double[maxNumDeps];
                    double deltaXi, deltaDFunc, deltaY, deltaYMag;
                    //double logDeltaXi, logDeltaD, logDeltaY, logDeltaYMag;

                    //Physical variables:
                    double[] radShell = new double[maxNumDeps]; //shell radial width                    
                    double[] tempShell = new double[maxNumDeps]; //kinetic temperature
                    double[] pressShell = new double[maxNumDeps]; //total pressure
                    double[] pGasShell = new double[maxNumDeps];  //gas pressure
                    double[] pRadShell = new double[maxNumDeps]; //radiation pressure
                    double[] rhoShell = new double[maxNumDeps];  //total fluid density
                    double[] massShell = new double[maxNumDeps];   // shell mass
                    double[] lumShell = new double[maxNumDeps]; //shell luminosity
                    double[] lumPpShell = new double[maxNumDeps]; //shell luminosity
                    double[] lumCnoShell = new double[maxNumDeps]; //shell luminosity                    

                    double[] epsShell = new double[maxNumDeps]; //total nuclear energy generation rate
                    double[] epsPpShell = new double[maxNumDeps]; //nuclear p-p chain energy generation rate  
                    double[] epsCnoShell = new double[maxNumDeps]; //nuclear cno cycle energy generation rate  

                    double[] kapShell = new double[maxNumDeps]; //mean opacity
                    double[] kapBfShell = new double[maxNumDeps]; //mean b-f opacity
                    double[] kapFfShell = new double[maxNumDeps]; //mean f-f opacity
                    double[] kapEsShell = new double[maxNumDeps]; //mean e^- scattering opacity
                    double[] kapHminShell = new double[maxNumDeps]; //mean H^- opacity                    

                    //cumulative quantities
                    double[] radInt = new double[maxNumDeps];  //interior radius
                    //double[] logRadInt = new double[maxNumDeps];
                    double[] massInt = new double[maxNumDeps];   //interior mass
                    double[] lumInt = new double[maxNumDeps]; //interior luminosity
                    double[] lumPpInt = new double[maxNumDeps]; //interior luminosity   
                    double[] lumCnoInt = new double[maxNumDeps]; //interior luminosity                                
                    double[] gravInt = new double[maxNumDeps]; //acceleration of gravity

                    //Convection:
                    double[] dLnPdLnT = new double[maxNumDeps];
                    boolean[] convFlag = new boolean[maxNumDeps];

                    double RHS, logRHS, logRHSMag; //, logLastXi;
                    double xiSquare; //useful

                    //Try uniform spacing for now...
                    // We know that at surface, xi >~ 3.0
                    //guess at a good spacing for now...
                    // Surface value of dimensionless xi parameter is in 
                    //   range 3 to 7 for polytropic index 1.5 to 3.0
                    deltaXi = 10.0 / maxNumDeps;
                    //deltaXi = 0.04;  //debug mode
                    //logDeltaXi = Math.log(deltaXi);

                    //For Newton-Raphson temperature recovery:
                    double firstTemp = 0.0;

                    //central bounday (initial) values:
                    // NOTE: we cannot set xi=0 - singularity                    
                    // We know that at surface, xi >~ 3.0
                    //guess at a good initial abscissa and spacing for now...
                    int j = 0;
                    xi[j] = 0.001;
                    yFunc[j] = 0.0;
                    DFunc[j] = 1.0;
                    // The stuff that follows...
                    radInt[j] = lambda * xi[j];
                    radShell[j] = radInt[j];
                    rhoShell[j] = rhoC * Math.pow(DFunc[j], index);
                    pressShell[j] = Kay * Math.pow(rhoShell[j], gammaPoly);
                    tempShell[j] = NRtemp.getTemp(pressShell[j], rhoShell[j], xFrac, yFrac, zFrac, firstTemp);
                    //pGasShell[j] = beta * pressShell[j];
                    //pRadShell[j] = (1.0 - beta) * pressShell[j];
                    pGasShell[j] = Useful.k * rhoShell[j] * tempShell[j] / PhysData.muI(xFrac, yFrac, zFrac) / Useful.amu;
                    pRadShell[j] = Useful.aStef() * Math.pow(tempShell[j], 4.0) / 3.0;
                    //tempShell[i] = pressShell[i] * PhysData.muI(xFrac, yFrac, zFrac) * Useful.amu / Useful.k / rhoShell[i];
                    massInt[j] = rhoC * 4.0 * Math.PI * Math.pow(radInt[j], 3) / 3.0;
                    massShell[j] = massInt[j];

                    if (tempShell[j] >= PhysData.fusionPPTemp) {
                        epsPpShell[j] = Power.ppChain(tempShell[j], rhoShell[j], xFrac, zFrac); //H fusion p-p chain
                        epsCnoShell[j] = Power.cnoCycle(tempShell[j], rhoShell[j], xFrac, zFrac); //H fusion CNO cycle
                        epsShell[j] = epsPpShell[j] + epsCnoShell[j];
                        //lumInt[j] = rhoC * 4.0 * Math.PI * Math.pow(radInt[j], 3) * epsShell[j] / 3.0;
                        lumPpInt[j] = rhoC * 4.0 * Math.PI * Math.pow(radInt[j], 3) * epsPpShell[j] / 3.0;
                        lumCnoInt[j] = rhoC * 4.0 * Math.PI * Math.pow(radInt[j], 3) * epsCnoShell[j] / 3.0;
                        lumPpShell[j] = lumPpInt[j];
                        lumCnoShell[j] = lumCnoInt[j];
                        lumInt[j] = lumPpInt[j] + lumCnoInt[j];
                        lumShell[j] = lumPpShell[j] + lumCnoShell[j];
                    } else {
                        epsPpShell[j] = 0.0; //H fusion p-p chain
                        epsCnoShell[j] = 0.0; //H fusion CNO cycle
                        epsShell[j] = 0.0;
                        //lumInt[j] = rhoC * 4.0 * Math.PI * Math.pow(radInt[j], 3) * epsShell[j] / 3.0;
                        lumPpInt[j] = 0.0;
                        lumCnoInt[j] = 0.0;
                        lumPpShell[j] = 0.0;
                        lumCnoShell[j] = 0.0;
                        lumInt[j] = 0.0;
                        lumShell[j] = 0.0;
                    }

                    kapBfShell[j] = Kappa.kappaBfFn(tempShell[j], rhoShell[j], xFrac, zFrac); //b-f photo-ionization
                    kapFfShell[j] = Kappa.kappaFfFn(tempShell[j], rhoShell[j], xFrac, zFrac); //f-f Bremsstrahlung
                    kapEsShell[j] = Kappa.kappaEsFn(tempShell[j], rhoShell[j], xFrac, zFrac); //Thomson e^- scattering
                    kapHminShell[j] = Kappa.kappaHminFn(tempShell[j], rhoShell[j], xFrac, zFrac); //H^- b-f 
                    kapShell[j] = kapBfShell[j] + kapFfShell[j] + kapEsShell[j] + kapHminShell[j];

                    gravInt[j] = Useful.GConst * massInt[j] / Math.pow(radInt[j], 2);

                    dLnPdLnT[j] = 0.0;

                    //4th order Runge-Kutta (RK4) helper variables
                    double k1y, k2y, k3y, k4y;
                    double k1D, k2D, k3D, k4D;
                    double hHalf = deltaXi / 2.0;
                    double yHalf, DHalf, xiHalf;
                    double yFull, DFull, xiFull;

                    //System.out.println("    i     " + "   radius   " + "   massIn    " + "   lumInt   " + "   temp   " + "   press   "
                    //        + "   rh    " + "   kappa   " + "   epsilon   " + "   dLnPdLnT    " + "   gravInt   " + "   pGas/PTot    "
                    //        + "   convection? ");
                    //System.out.println("      i       " + "    radius    " + "    press    " + "    rho    ");
                    //Master numerical integration loop:
                    //Loop exits when surface is found
                    // Try Euler's method for now...
                    //System.out.println("i    xi[i]    yFunc[i]    DFunc[i] ");
                    int iSurf = 0;
                    int iCore = 0;
                    for (int i = 1; i < maxNumDeps; i++) {

                        xi[i] = xi[i - 1] + deltaXi;

                        /* Logarithmic - difficult due to signs
                         logLastXi = Math.log(xi[i - 1]);

                         //Start by advanceing helper function, y(xi):
                         logRHSMag = index * Math.log(DFunc[i - 1]) + 2.0 * logLastXi;
                         //RHS = -1.0 * Math.exp(logRHSMag);
                         logDeltaYMag = logRHSMag + logDeltaXi;
                         deltaYMag = Math.exp(logDeltaYMag);
                         deltaY = -1.0 * deltaYMag;
                         yFunc[i] = yFunc[i - 1] + deltaY;

                         //Now use the updated value of y to update D_n
                         logRHS = Math.log(yFunc[i]) - 2.0 * logLastXi;
                         logDeltaD = logRHS + logDeltaXi;
                         deltaDFunc = Math.exp(logDeltaD);
                         DFunc[i] = DFunc[i - 1] + deltaDFunc;
                         */
                        //
                        xiSquare = Math.pow(xi[i - 1], 2);
                        xiHalf = Math.pow((xi[i - 1] + hHalf), 2);
                        xiFull = Math.pow(xi[i], 2);
                        //

                        /*
                         //Euler's method
                         //Start by advanceing helper function, y(xi):
                         RHS = -1.0 * Math.pow(DFunc[i - 1], index) * xiSquare;
                         deltaY = RHS * deltaXi;
                         yFunc[i] = yFunc[i - 1] + deltaY;

                         //Now use the updated value of y to update D_n
                         RHS = yFunc[i] / xiSquare;
                         deltaDFunc = RHS * deltaXi;
                         DFunc[i] = DFunc[i - 1] + deltaDFunc;
                        
                         System.out.format("Euler: %03d   %15.10f   %15.10f   %15.10f%n", i, xi[i], yFunc[i], DFunc[i]);
                         */
                        //4th order Runge-Kutta (RK4):
                        k1y = -1.0 * Math.pow(DFunc[i - 1], index) * xiSquare;
                        k1D = yFunc[i - 1] / xiSquare;

                        DHalf = DFunc[i - 1] + hHalf * k1D;
                        yHalf = yFunc[i - 1] + hHalf * k1y;
                        k2y = -1.0 * Math.pow(DHalf, index) * xiHalf;
                        k2D = yHalf / xiHalf;

                        DHalf = DFunc[i - 1] + hHalf * k2D;
                        yHalf = yFunc[i - 1] + hHalf * k2y;
                        k3y = -1.0 * Math.pow(DHalf, index) * xiHalf;
                        k3D = yHalf / xiHalf;

                        DFull = DFunc[i - 1] + deltaXi * k3D;
                        yFull = yFunc[i - 1] + deltaXi * k3y;
                        k4y = -1.0 * Math.pow(DFull, index) * xiFull;
                        k4D = yFull / xiFull;

                        deltaY = deltaXi * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
                        deltaDFunc = deltaXi * (k1D + 2.0 * k2D + 2.0 * k3D + k4D) / 6.0;
                        yFunc[i] = yFunc[i - 1] + deltaY;
                        DFunc[i] = DFunc[i - 1] + deltaDFunc;

                        //Are we there yet?
                        if ((DFunc[i] <= 0)
                                || (Double.isNaN(DFunc[i]) == true)
                                || (Double.isInfinite(DFunc[i]) == true)) {
                            break;
                        }

                        //System.out.format("RK4: %03d   %15.10f   %15.10f   %15.10f%n", i, xi[i], yFunc[i], DFunc[i]);
                        radInt[i] = lambda * xi[i];
                        radShell[i] = radInt[i] - radInt[i - 1];
                        rhoShell[i] = rhoC * Math.pow(DFunc[i], index);
                        pressShell[i] = Kay * Math.pow(rhoShell[i], gammaPoly);
                        tempShell[i] = NRtemp.getTemp(pressShell[i], rhoShell[i], xFrac, yFrac, zFrac, tempShell[i - 1]);
                        //pGasShell[i] = beta * pressShell[i];
                        //pRadShell[i] = (1.0 - beta) * pressShell[i];
                        pGasShell[i] = Useful.k * rhoShell[i] * tempShell[i] / PhysData.muI(xFrac, yFrac, zFrac) / Useful.amu;
                        pRadShell[i] = Useful.aStef() * Math.pow(tempShell[i], 4.0) / 3.0;
                        //tempShell[i] = pressShell[i] * PhysData.muI(xFrac, yFrac, zFrac) * Useful.amu / Useful.k / rhoShell[i];
                        massShell[i] = 4.0 * Math.PI * Math.pow(radInt[i], 2) * rhoShell[i] * radShell[i];
                        massInt[i] = massInt[i - 1] + massShell[i];
                        //epsShell[i] = Power.nuclear(tempShell[i], rhoShell[i], xFrac, zFrac);

                        if (tempShell[i] >= PhysData.fusionPPTemp) {
                            iCore = i;
                            epsPpShell[i] = Power.ppChain(tempShell[i], rhoShell[i], xFrac, zFrac); //H fusion p-p chain
                            epsCnoShell[i] = Power.cnoCycle(tempShell[i], rhoShell[i], xFrac, zFrac); //H fusion CNO cycle
                            epsShell[i] = epsPpShell[i] + epsCnoShell[i];
                            //lumShell[i] = 4.0 * Math.PI * Math.pow(radInt[i], 2) * rhoShell[i] * epsShell[i] * radShell[i];
                            lumPpShell[i] = 4.0 * Math.PI * Math.pow(radInt[i], 2) * rhoShell[i] * epsPpShell[i] * radShell[i];
                            lumCnoShell[i] = 4.0 * Math.PI * Math.pow(radInt[i], 2) * rhoShell[i] * epsCnoShell[i] * radShell[i];
                            lumShell[i] = lumPpShell[i] + lumCnoShell[i];
                            lumPpInt[i] = lumPpInt[i - 1] + lumPpShell[i];
                            lumCnoInt[i] = lumCnoInt[i - 1] + lumCnoShell[i];
                            lumInt[i] = lumInt[i - 1] + lumShell[i];
                        } else {
                            epsPpShell[i] = 0.0; //H fusion p-p chain
                            epsCnoShell[i] = 0.0; //H fusion CNO cycle
                            epsShell[i] = 0.0;
                            //lumShell[i] = 4.0 * Math.PI * Math.pow(radInt[i], 2) * rhoShell[i] * epsShell[i] * radShell[i];
                            lumPpShell[i] = 0.0;
                            lumCnoShell[i] = 4.0 * Math.PI * Math.pow(radInt[i], 2) * rhoShell[i] * epsCnoShell[i] * radShell[i];
                            lumShell[i] = 0.0;
                            lumPpInt[i] = lumPpInt[i - 1];
                            lumCnoInt[i] = lumCnoInt[i - 1];
                            lumInt[i] = lumInt[i - 1];
                        }

                        //kapShell[i] = Kappa.kappaFn(tempShell[i], rhoShell[i], xFrac, zFrac);
                        kapBfShell[i] = Kappa.kappaBfFn(tempShell[i], rhoShell[i], xFrac, zFrac); //b-f photo-ionization
                        kapFfShell[i] = Kappa.kappaFfFn(tempShell[i], rhoShell[i], xFrac, zFrac); //f-f Bremsstrahlung
                        kapEsShell[i] = Kappa.kappaEsFn(tempShell[i], rhoShell[i], xFrac, zFrac); //Thomson e^- scattering
                        kapHminShell[i] = Kappa.kappaHminFn(tempShell[i], rhoShell[i], xFrac, zFrac); //H^- b-f
                        kapShell[i] = kapBfShell[i] + kapFfShell[i] + kapEsShell[i] + kapHminShell[i];

                        gravInt[i] = Useful.GConst * massInt[i] / Math.pow(radInt[i], 2);

                        dLnPdLnT[i] = (Math.log(pGasShell[i]) - Math.log(pGasShell[i - 1]))
                                / (Math.log(tempShell[i]) - Math.log(tempShell[i - 1]));
                        if ((dLnPdLnT[i] >= gamThing)) {
                            //Radiative transport
                            convFlag[i] = false;
                        } else {
                            //Convective transport
                            convFlag[i] = true;
                        }

                        //Dynamically update Kay and lambda:
                        //Hmmm... these never seem to change, but let's update them anyway...
                        Kay = pressShell[i] / Math.pow(rhoShell[i], gammaPoly);
                        lamTerm = (index + 1.0) * Kay * Math.pow(rhoC, indxExp) / 4.0 / Math.PI / Useful.GConst;
                        lambda = Math.pow(lamTerm, 0.5);
                        //System.out.println("Kay " + Kay + " lambda " + lambda);

                        //System.out.println("    i     " + "     radius     " + "     massInt     " + "     lumInt     " + "     temp     " + "     press     "
                        //         + "     rho     " + "     kappa     " + "     epsilon     " + "     dLnPdLnT     " + "     gravInt     " + "     pGas/PTot     " 
                        //         + "     convection? ");
                        //System.out.println("      i       " + "    radius    " + "    press    " + "    rho    ");                        
                        //System.out.format("%03d  %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %b%n",
                        //        i, log10E * Math.log(radInt[i]), log10E * Math.log(massInt[i]), log10E * Math.log(lumInt[i]), log10E * Math.log(tempShell[i]),
                        //        log10E * Math.log(pressShell[i]), log10E * Math.log(rhoShell[i]), kapShell[i], log10E * Math.log(epsShell[i]),
                        //        dLnPdLnT[i], log10E * Math.log(gravInt[i]), pGasShell[i] / pressShell[i], convFlag[i]);
                        //System.out.format("   %03d     %15.11f      %15.11f      %15.11f      %15.11f      %n",
                        //        i, log10E * Math.log(radInt[i]), log10E * Math.log(pressShell[i]), log10E * Math.log(rhoShell[i]), log10E * Math.log(tempShell[i]));
                        //Surface boundary condition:
                        iSurf++;

                    }

                    int numDeps = iSurf;
                    System.out.println("Actual number of depths = " + numDeps);
                    //Independent total mass calculation:
                    //First derivative, dD_n/dXi at surface:
                    double dDFuncdXi = (DFunc[iSurf] - DFunc[iSurf - 1]) / (xi[iSurf] - xi[iSurf - 1]);
                    double totMass = -4.0 * Math.PI * Math.pow(lambda, 3.0) * rhoC * Math.pow(xi[iSurf], 2.0) * dDFuncdXi;
                    System.out.println("Total analytic mass " + totMass);

                    //Results:
                    //cgs units:
                    double mass = massInt[numDeps];
                    double radius = radInt[numDeps];
                    double luminosity = lumInt[numDeps];
                    double surfTemp = tempShell[numDeps];

                    //Compute the effective temperature of the model:
                    double teff4 = luminosity / 4.0 / Math.PI / Math.pow(radius, 2.0) / Useful.sigma;
                    double teff = Math.pow(teff4, 0.25);

                    //solar units:
                    double massSol = mass / Useful.mSun;
                    double radiusSol = radius / Useful.rSun;
                    double lumSol = luminosity / Useful.lSun;

                    //Run the loop inward to build the optical depth scasle, tauIn:
                    //// No!While we're at it - find the nuclear burning core
                    double[] tauIn = new double[numDeps];
                    tauIn[numDeps - 1] = 0.0;
                    //int iCore = numDeps - 1;
                    for (int i = numDeps - 2; i == 0; i--) {
                        tauIn[i] = tauIn[i - 1]
                                + rhoShell[i] * kapShell[i] * radShell[i];
                       // if (lumInt[i] > 0.99 * luminosity) {
                        //    iCore--;
                       // }
                    }

                    //Report:
                    System.out.println("cgs units:");
                    System.out.println("Radius: " + log10E * Math.log(radius) + " Mass: " + log10E * Math.log(mass) + " Bol Luminosity: " + log10E * Math.log(luminosity));
                    System.out.println("Teff: " + teff + " TSurf: " + surfTemp);
                    System.out.println("Solar units:");
                    System.out.println("Radius: " + radiusSol + " Mass: " + massSol + " Bol Luminosity: " + lumSol);
                    System.out.println("Nuclear burning core fractional radius: " + (radInt[iCore] / radius));

                    //Compute spectral energy distribution (SED):
                    //wavelength grid (cm):
                    double[] waveSetup = new double[3];
                    waveSetup[0] = 100.0 * 1.0e-7;  // test Start wavelength, cm
                    waveSetup[1] = 2000.0 * 1.0e-7; // test End wavelength, cm
                    waveSetup[2] = 100;  // test number of lambda
                    int numWaves = (int) waveSetup[2];

                    double[] SED = new double[numWaves];
                    double[] waveGrid = new double[numWaves];
                    double thisWave, thisLogWave, logWave0, logWave1;
                    logWave0 = Math.log(waveSetup[0]);
                    logWave1 = Math.log(waveSetup[1]);
                    double deltaLogWave = (logWave1 - logWave0) / numWaves;
                    for (int i = 0; i < numWaves; i++) {
                        thisLogWave = logWave0 + ((double) i) * deltaLogWave;
                        thisWave = Math.exp(thisLogWave);
                        waveGrid[i] = thisWave;
                        SED[i] = Planck.planck(teff, thisWave);
                        //System.out.println(" " + waveGrid[i] + " " + SED[i]);
                    }

                    //
                    //double colors[] = new double[5];
                    double colors[] = Photometry.UBVRI(waveGrid, SED);

                    // All code before this line
                    String patternCol = "0.00";
                    //String pattern = "#####.##";
                    DecimalFormat colFormatter = new DecimalFormat(patternCol);
                    // // String patternWl = "0.00";
                    //  //String pattern = "#####.##";
                    // DecimalFormat WlFormatter = new DecimalFormat(patternWl);

                    actiontarget.setText("Photometric color indices: "
                            + "U-B: " + colFormatter.format(colors[0])
                            + " B-V: " + colFormatter.format(colors[1])
                            + " V-R: " + colFormatter.format(colors[2])
                            + " V-I: " + colFormatter.format(colors[3])
                            + " R-I: " + colFormatter.format(colors[4]) + "\r\n"
                    );

                    // No! //grid.getChildren().add(r);
                    // Graphical output section:
                    /*
                    // Plot 1: T_Kin(radius):
                    LineChart<Number, Number> lineChartT2 = LineCharts.t2Plot(numDeps, radInt, tempShell);
                    grid.add(lineChartT2, 0, 9);
                    // Plot 2: log(P(radius):
                    LineChart<Number, Number> lineChartP = LineCharts.pressPlot(numDeps, radInt, pressShell, pGasShell, pRadShell);
                    grid.add(lineChartP, 2, 9);
                    // Plot 3: log rho:
                    LineChart<Number, Number> lineChartRho = LineCharts.rhoPlot(numDeps, radInt, rhoShell);
                    grid.add(lineChartRho, 2, 11);

                    // Plot 4: log cumulative L_Bol
                    LineChart<Number, Number> lineChartLum = LineCharts.lumPlot(numDeps, radInt, lumInt, lumPpInt, lumCnoInt);
                    grid.add(lineChartLum, 0, 10);

                    // Plot 5: log cumulative mass
                    LineChart<Number, Number> lineChartMass = LineCharts.massPlot(numDeps, radInt, massInt);
                    grid.add(lineChartMass, 2, 10);
                    // Plot 6: log epsilon (nuc power generation)
                    LineChart<Number, Number> lineChartEps = LineCharts.epsPlot(numDeps, radInt, epsShell, epsPpShell, epsCnoShell);
                    grid.add(lineChartEps, 0, 11);
                    // Plot 7: log kappa (mass extinction)
                    LineChart<Number, Number> lineChartKap = LineCharts.kapPlot(numDeps, radInt, kapShell, kapBfShell, kapFfShell, kapEsShell, kapHminShell);
                    grid.add(lineChartKap, 4, 10);
                    // Plot 8: log flux(lambda) (SED)
                    LineChart<Number, Number> lineChartSED = LineCharts.sedPlot(numWaves, waveGrid, SED);
                    grid.add(lineChartSED, 4, 11);
                    */
                    
                    // Graphical output section:
                    // Plot 1: T_Kin(radius):
                    LineChart<Number, Number> lineChartT2 = LineCharts.t2Plot(numDeps, radInt, tempShell);
                    grid.add(lineChartT2, 0, 9);
                    // Plot 2: log(P(radius):
                    LineChart<Number, Number> lineChartP = LineCharts.pressPlot(numDeps, radInt, pressShell, pGasShell, pRadShell);
                    grid.add(lineChartP, 2, 9);
                    // Plot 3: log rho:
                    LineChart<Number, Number> lineChartRho = LineCharts.rhoPlot(numDeps, radInt, rhoShell);
                    grid.add(lineChartRho, 2, 13);

                    // Plot 4: log cumulative L_Bol
                    LineChart<Number, Number> lineChartLum = LineCharts.lumPlot(numDeps, iCore, radInt, lumInt, lumPpInt, lumCnoInt);
                    grid.add(lineChartLum, 0, 11);

                    // Plot 5: log cumulative mass
                    LineChart<Number, Number> lineChartMass = LineCharts.massPlot(numDeps, radInt, massInt);
                    grid.add(lineChartMass, 2, 11);
                    // Plot 6: log epsilon (nuc power generation)
                    LineChart<Number, Number> lineChartEps = LineCharts.epsPlot(numDeps, radInt, epsShell, epsPpShell, epsCnoShell);
                    grid.add(lineChartEps, 0, 13);
                    // Plot 7: log kappa (mass extinction)
                    LineChart<Number, Number> lineChartKap = LineCharts.kapPlot(numDeps, radInt, kapShell, kapBfShell, kapFfShell, kapEsShell, kapHminShell);
                    grid.add(lineChartKap, 4, 11);
                    // Plot 8: log flux(lambda) (SED)
                    LineChart<Number, Number> lineChartSED = LineCharts.sedPlot(numWaves, waveGrid, SED);
                    grid.add(lineChartSED, 4, 13);
                    

                    //                   LineChart<Number, Number> lineChartSpec = LineCharts.specPlot(numMaster, masterLams, cosTheta, masterIntens, masterFlux);
                    //                   grid.add(lineChartSpec, 2, 10);
                    ////LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, lineProf, lam0);
                    //                   grid.add(lineChartLine, 4, 10);
                    //Debug versions of the function with parameters for plotting up quqntities used in line profile
                    // calculation:
                    //LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, logKappaL, lam0, kappa);
                    //grid.add(lineChartLine, 4, 10);
                    //LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, logTauL, lam0);
                    //grid.add(lineChartLine, 4, 10);
                    // Final scene / stage stuff:
                    //Scene scene = new Scene(grid, 1600, 900);
                    //scene.getStylesheets().add("../../GrayCascadeStyleSheet.css");
                    //
                    //primaryStage.setScene(scene);
                    //
                    //primaryStage.show();
                } else {
                    actiontarget.setText("All fields must have values");
                }

            }
        }
        );

        // Final scene / stage stuff:
        Scene scene = new Scene(grid, 1600, 900);

        //Stylesheet.css must go in /src/ directory
        scene.getStylesheets()
                .add("GrayCascadeStyleSheet.css");

        primaryStage.setScene(scene);
        primaryStage.show();
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        launch(args);
    }

}
