/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package polystar;

import java.text.DecimalFormat;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
//import javafx.scene.chart.XYChart.Data;
//import javafx.scene.chart.XYChart.getXValue;

/**
 *
 * @author Ian
 */
public class LineCharts {

    //Plot one: T_Kin vs radius
    public static LineChart<Number, Number> t2Plot(int numDeps, double[] depths, double[] temp) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((0.9 * temp[0]) / 1000.0);
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((1.1 * temp[numDeps - 1]) / 1000.0);
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "T_Kin (K)";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartT2;
        lineChartT2 = new LineChart<>(xAxis, yAxis);

        lineChartT2.setTitle("Kinetic temperature vs radius"); //\r\n" + teffLbl + "K");
        lineChartT2.setId("tempDpth");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesT;
        seriesT = new XYChart.Series();
        seriesT.setName("T_Kin");

        double depths2[] = new double[numDeps];
        double conv = 1.0E-5;
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            XYChart.Data<?, ?> xyPair;
            xyPair = new XYChart.Data(depths2[id], temp[id]);
            boolean add;
            //add = seriesT.getData().add(new XYChart.Data(conv * depths[id], temp[0][id]));
            add = seriesT.getData().add(xyPair);
        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartT2.getData().addAll(seriesT); //, seriesdTau1, seriestTau1);

        lineChartT2.setCreateSymbols(false);

        return lineChartT2;

    }

    //Plot two: Press vs radius with pGas and pRad
    public static LineChart<Number, Number> pressPlot(int numDeps, double[] depths, double[] press, double[] pGas, double[] pRad) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((logE * Math.log(0.9 * press[0]) / 1000.0));
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((logE * Math.log(1.1 * press[numDeps - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 Pressure";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartP;
        lineChartP = new LineChart<>(xAxis, yAxis);

        lineChartP.setTitle("log_10 Pressure vs radius"); //\r\n" + teffLbl + "K");
        lineChartP.setId("pressDpth");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesP = new XYChart.Series();
        seriesP.setName("logP_Tot");

        XYChart.Series seriesPg = new XYChart.Series();
        seriesPg.setName("logP_Gas");

        XYChart.Series seriesPr = new XYChart.Series();
        seriesPr.setName("logP_Rad");

        double[] depths2 = new double[numDeps];
        double conv = 1.0E-5;

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            seriesP.getData().add(new XYChart.Data(depths2[id], logE * Math.log(press[id]))); // Total pressure
            seriesPg.getData().add(new XYChart.Data(depths2[id], logE * Math.log(pGas[id]))); // gas pressure
            seriesPr.getData().add(new XYChart.Data(depths2[id], logE * Math.log(pRad[id]))); // radiation pressure

        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartP.getData().addAll(seriesP, seriesPg, seriesPr); //, seriesdTau1, seriestTau1);

        lineChartP.setCreateSymbols(false);

        return lineChartP;

    }

    //Plot three: rho vs radius 
    public static LineChart<Number, Number> rhoPlot(int numDeps, double[] depths, double[] rho) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((logE * Math.log(0.9 * rho[0]) / 1000.0));
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((logE * Math.log(1.1 * rho[numDeps - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 rho";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartRho;
        lineChartRho = new LineChart<>(xAxis, yAxis);

        lineChartRho.setTitle("log_10 mass density vs radius");
        lineChartRho.setId("logRho");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesRho = new XYChart.Series();
        seriesRho.setName("logRho");

        double[] depths2 = new double[numDeps];
        double conv = 1.0E-5;

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            seriesRho.getData().add(new XYChart.Data(depths2[id], logE * Math.log(rho[id]))); // Mass density

        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartRho.getData().addAll(seriesRho); //, seriesdTau1, seriestTau1);

        lineChartRho.setCreateSymbols(false);

        return lineChartRho;

    }

    //Plot four: cumulative luminosity vs radius 
    public static LineChart<Number, Number> lumPlot(int numDeps, int iCore, double[] depths, double[] lumInt, double[] lumPpInt, double[] lumCnoInt) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((logE * Math.log(0.9 * lumInt[0]) / 1000.0));
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((logE * Math.log(1.1 * lumInt[numDeps - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;
        
        double conv = 1.0E-5;        
        
        //Nuclear burning core boundary marker:
        double rCore = conv * depths[iCore];
        //double[] xCore = new double[2];
        double[] xCore = {rCore, rCore};
        double[] yCore = {minY, 40+maxY}; // I don't know what's happening with maxY!

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 L_Bol";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartLum;
        lineChartLum = new LineChart<>(xAxis, yAxis);

        lineChartLum.setTitle("log_10 Cumulative Bolometric luminosity vs radius");
        lineChartLum.setId("logLum");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesLum = new XYChart.Series();
        seriesLum.setName("logLum");
        XYChart.Series seriesLumPp = new XYChart.Series();
        seriesLumPp.setName("logLumPp");
        XYChart.Series seriesLumCno = new XYChart.Series();
        seriesLumCno.setName("logLumCno");
        XYChart.Series seriesCore = new XYChart.Series();
        seriesCore.setName("core");

        double[] depths2 = new double[numDeps];


        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            seriesLum.getData().add(new XYChart.Data(depths2[id], logE * Math.log(lumInt[id]))); // total luminosity
            seriesLumPp.getData().add(new XYChart.Data(depths2[id], logE * Math.log(lumPpInt[id]))); // p-p chain luminosity
            seriesLumCno.getData().add(new XYChart.Data(depths2[id], logE * Math.log(lumCnoInt[id]))); // cno chain luminosity

        }
        
        //nuclear burning core boundary:
        seriesCore.getData().add(new XYChart.Data(xCore[0], yCore[0]));
        seriesCore.getData().add(new XYChart.Data(xCore[1], yCore[1]));        

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartLum.getData().addAll(seriesLum, seriesLumPp, seriesLumCno, seriesCore); //, seriesdTau1, seriestTau1);

        lineChartLum.setCreateSymbols(false);

        return lineChartLum;

    }

    //Plot five: cumulative mass vs radius 
    public static LineChart<Number, Number> massPlot(int numDeps, double[] depths, double[] massInt) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((logE * Math.log(0.9 * massInt[0]) / 1000.0));
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((logE * Math.log(1.1 * massInt[numDeps - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 Mass";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartMass;
        lineChartMass = new LineChart<>(xAxis, yAxis);

        lineChartMass.setTitle("log_10 Cumulative mass vs radius"); //\r\n" + teffLbl + "K");
        lineChartMass.setId("logMass");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesMass = new XYChart.Series();
        seriesMass.setName("logMass");

        double[] depths2 = new double[numDeps];
        double conv = 1.0E-5;

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            seriesMass.getData().add(new XYChart.Data(depths2[id], logE * Math.log(massInt[id]))); // interior mass 

        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartMass.getData().addAll(seriesMass); //, seriesdTau1, seriestTau1);

        lineChartMass.setCreateSymbols(false);

        return lineChartMass;

    }

    //Plot six: epsilon - nuclear power generation 
    public static LineChart<Number, Number> epsPlot(int numDeps, double[] depths, double[] epsShell, double[] epsPpShell, double[] epsCnoShell) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((logE * Math.log(0.9 * epsShell[0]) / 1000.0));
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((logE * Math.log(1.1 * epsShell[numDeps - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 epsilon";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartEps;
        lineChartEps = new LineChart<>(xAxis, yAxis);

        lineChartEps.setTitle("log_10 epsilon vs radius"); //\r\n" + teffLbl + "K");
        lineChartEps.setId("logEps");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesEps = new XYChart.Series();
        seriesEps.setName("logEps");
        XYChart.Series seriesEpsPp = new XYChart.Series();
        seriesEpsPp.setName("logEpsPp");
        XYChart.Series seriesEpsCno = new XYChart.Series();
        seriesEpsCno.setName("logEpsCno");

        double[] depths2 = new double[numDeps];
        double conv = 1.0E-5;

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            if (epsShell[id] <= 0.0) {
                seriesEps.getData().add(new XYChart.Data(depths2[id], -9.0)); // total H fusion power generation
                seriesEpsPp.getData().add(new XYChart.Data(depths2[id], -9.0)); // p-p chain power generation
                seriesEpsCno.getData().add(new XYChart.Data(depths2[id], -9.0)); // CNO cycle power generation
            } else {
                seriesEps.getData().add(new XYChart.Data(depths2[id], logE * Math.log(epsShell[id]))); // total H fusion power generation
                seriesEpsPp.getData().add(new XYChart.Data(depths2[id], logE * Math.log(epsPpShell[id]))); // p-p chain power generation
                seriesEpsCno.getData().add(new XYChart.Data(depths2[id], logE * Math.log(epsCnoShell[id]))); // CNO cycle power generation                
            }

        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartEps.getData().addAll(seriesEps, seriesEpsPp, seriesEpsCno);

        lineChartEps.setCreateSymbols(false);

        return lineChartEps;

    }

    //Plot seven: kappa - opacity 
    public static LineChart<Number, Number> kapPlot(int numDeps, double[] depths, double[] kapShell,
            double[] kapBfShell, double[] kapFfShell, double[] kapEsShell, double[] kapHminShell) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((logE * Math.log(0.9 * kapShell[0]) / 1000.0));
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((logE * Math.log(1.1 * kapShell[numDeps - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Radius (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 kappa";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartKap;
        lineChartKap = new LineChart<>(xAxis, yAxis);

        lineChartKap.setTitle("log_10 kappa vs radius"); //\r\n" + teffLbl + "K");
        lineChartKap.setId("logKap");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesKap = new XYChart.Series();
        seriesKap.setName("logKap");
        XYChart.Series seriesKapBf = new XYChart.Series();
        seriesKapBf.setName("logKapBf");
        XYChart.Series seriesKapFf = new XYChart.Series();
        seriesKapFf.setName("logKapFf");
        XYChart.Series seriesKapEs = new XYChart.Series();
        seriesKapEs.setName("logKapEs");
        XYChart.Series seriesKapHmin = new XYChart.Series();
        seriesKapHmin.setName("logKapHmin");

        double[] depths2 = new double[numDeps];
        double conv = 1.0E-5;

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            seriesKap.getData().add(new XYChart.Data(depths2[id], logE * Math.log(kapShell[id]))); // mean mass extinction
            seriesKapBf.getData().add(new XYChart.Data(depths2[id], logE * Math.log(kapBfShell[id]))); // b-f contribution
            seriesKapFf.getData().add(new XYChart.Data(depths2[id], logE * Math.log(kapFfShell[id]))); // f-f (Bremsstrahlung) contribution
            seriesKapEs.getData().add(new XYChart.Data(depths2[id], logE * Math.log(kapEsShell[id]))); // Thomson e^- scattering contribution
            //seriesKapHmin.getData().add(new XYChart.Data(depths2[id], logE*Math.log(kapHminShell[id]))); // H^- b-f contribution  
            if (kapHminShell[id] <= 0) {
                seriesKapHmin.getData().add(new XYChart.Data(depths2[id], -9.0));
            } else {
                seriesKapHmin.getData().add(new XYChart.Data(depths2[id], logE * Math.log(kapHminShell[id]))); // H^- b-f contribution 
            }

        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartKap.getData().addAll(seriesKap, seriesKapBf, seriesKapFf, seriesKapEs, seriesKapHmin); //, seriesdTau1, seriestTau1);

        lineChartKap.setCreateSymbols(false);

        return lineChartKap;

    }

    //Plot five: cumulative mass vs radius 
    public static LineChart<Number, Number> sedPlot(int numWaves, double[] waves, double[] flux) {

        double logE = Math.log10(Math.E);

        double minX = waves[0];
        double maxX = waves[numWaves - 1];
        //double deltaX = 100.0;

        int minYint = (int) ((0.9 * flux[0]) / 1000.0);
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) (((1.1 * flux[numWaves - 1]) / 1000.0));
        double maxY = (double) (maxYint) * 1000.0;
//        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "Wavelength (nm)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "log_10 Flux";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartSED;
        lineChartSED = new LineChart<>(xAxis, yAxis);

        lineChartSED.setTitle("log_10 Flux vs lambda"); //\r\n" + teffLbl + "K");
        lineChartSED.setId("logSED");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesSED = new XYChart.Series();
        seriesSED.setName("logFlux");

        double[] waves2 = new double[numWaves];
        double conv = 1.0E7;
        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numWaves; id++) {

            waves2[id] = conv * waves[id];
            seriesSED.getData().add(new XYChart.Data(waves2[id], (flux[id]))); // SED

        }

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartSED.getData().addAll(seriesSED); //, seriesdTau1, seriestTau1);

        lineChartSED.setCreateSymbols(false);

        return lineChartSED;

    }

    /**
     * Plot 5: SED: log(F_lambda(log(lambda)), I_lambda(log(lambda))):
     *
     * @param numLams
     * @param lambdaScale
     * @param cosTheta
     * @param intens
     * @param flux
     * @param lamUBVRI
     * @return
     */
    /*
     public static LineChart<Number, Number> specPlot(int numLams,
     double[] lambdaScale, double[][] cosTheta, double[][] intens,
     double[][] flux) {

     double logE = Math.log10(Math.E);
     int numMus = cosTheta[0].length;

     final NumberAxis xAxis = new NumberAxis();
     final NumberAxis yAxis = new NumberAxis();

     xAxis.setLabel("lambda (nm)");

     String yTitle = "F_lam, I_lam (/F_max)";
     yAxis.setLabel(yTitle);

     final LineChart<Number, Number> lineChartSpec
     = new LineChart<Number, Number>(xAxis, yAxis);

     lineChartSpec.setId("spec");

     String pattern = "0.00";
     //String pattern = "###.####";
     DecimalFormat myFormatter = new DecimalFormat(pattern);

     XYChart.Series seriesF = new XYChart.Series();
     seriesF.setName("Flux");

     XYChart.Series seriesI0 = new XYChart.Series();
     seriesI0.setName("I(theta=0)");
     String thetaStr = myFormatter.format((180.0 / Math.PI) * Math.acos(cosTheta[1][numMus - 2]));
     XYChart.Series seriesIn = new XYChart.Series();
     seriesIn.setName("I(theta=" + thetaStr + ")");

     // Only linear units will work (!) - divide flux by lambda_max:
     int[] iLamMinMax = MinMax2.minMax(flux);
     int iLamMax = iLamMinMax[1];
     double fMax = flux[0][iLamMax];
     String lamStr = myFormatter.format(1.0e7 * lambdaScale[iLamMax]);

     lineChartSpec.setTitle("Spectral energy distribution (SED)\r\n" + "lambda_Max: " + lamStr + " nm");

     for (int il = 0; il < numLams - 1; il++) {

     lambdaScale[il] = lambdaScale[il] * 1.0e7;
     //System.out.println("il, lambdaScale[il], flux[0][il], logE*flux[1][il]: " + il + " " + lambdaScale[il] + " " + flux[0][il] + " " + logE * flux[1][il]);
     //seriesF.getData().add( new XYChart.Data( Math.log10(lambdaScale[il]), logE*flux[1][il] ) );
     seriesF.getData().add(new XYChart.Data(lambdaScale[il], flux[0][il] / fMax));
     //seriesI0.getData().add( new XYChart.Data( Math.log10(lambdaScale[il]), Math.log10(intens[il][0]) ) );
     //seriesIn.getData().add( new XYChart.Data( Math.log10(lambdaScale[il]), Math.log10(intens[il][numMus-2]) ) );
     seriesI0.getData().add(new XYChart.Data(lambdaScale[il], intens[il][0] / fMax));
     seriesIn.getData().add(new XYChart.Data(lambdaScale[il], intens[il][numMus - 2] / fMax));
     }

     //Add UBVRI band centres
     //Photometric bands - Johnson UBVRI
     double[][][] filters = FilterSet.filterSet();
     int lam0_ptr = 11; // approximate band centre
     int numBands = filters.length;
     double[] lamUBVRI = new double[numBands];
     for (int ib = 0; ib < numBands; ib++) {
     lamUBVRI[ib] = filters[ib][0][lam0_ptr] * 1.0e7;
     }
     //System.out.println("lamUBVRI[0]: " + lamUBVRI[0][1]);

     XYChart.Series seriesU = new XYChart.Series();
     seriesU.setName("U");
     seriesU.getData().add(new XYChart.Data(lamUBVRI[0], 0.0));
     seriesU.getData().add(new XYChart.Data(lamUBVRI[0], 1.1));
     XYChart.Series seriesB = new XYChart.Series();
     seriesB.setName("B");
     seriesB.getData().add(new XYChart.Data(lamUBVRI[2], 0.0));
     seriesB.getData().add(new XYChart.Data(lamUBVRI[2], 1.1));
     XYChart.Series seriesV = new XYChart.Series();
     seriesV.setName("V");
     seriesV.getData().add(new XYChart.Data(lamUBVRI[3], 0.0));
     seriesV.getData().add(new XYChart.Data(lamUBVRI[3], 1.1));
     XYChart.Series seriesR = new XYChart.Series();
     seriesR.setName("R");
     seriesR.getData().add(new XYChart.Data(lamUBVRI[4], 0.0));
     seriesR.getData().add(new XYChart.Data(lamUBVRI[4], 1.1));
     XYChart.Series seriesI = new XYChart.Series();
     seriesI.setName("I");
     seriesI.getData().add(new XYChart.Data(lamUBVRI[5], 0.0));
     seriesI.getData().add(new XYChart.Data(lamUBVRI[5], 1.1));

     //lineChartSpec.getData().add(seriesSpec);
     lineChartSpec.getData().addAll(seriesF, seriesI0, seriesIn,
     seriesU, seriesB, seriesV, seriesR, seriesI);

     lineChartSpec.setCreateSymbols(false);

     return lineChartSpec;

     } */
}
