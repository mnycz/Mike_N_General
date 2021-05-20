ackage org.clas.modules;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedList;

import java.util.ArrayList;
import java.util.List;

public class ALERTDataStructs implements IDataEventListener{
    public DataGroup SCDG;
    public static IndexedList<DataGroup> dataGroups = new IndexedList<DataGroup>(3);
    public static H1F Pedestal[][][] = new H1F [16][3][5];
    public static H1F fbAlign[][][] = new H1F [16][3][5];
    public static H2F TW[][][] = new H2F [16][3][5];
    public static H2F Veff[][][] = new H2F [16][3][5];

    //@Override


    public void dataEventAction(DataEvent dataEvent) {
        System.out.println("In Data Struct");
        ALERTCalibrationEngine CalibrationRoutines = new ALERTCalibrationEngine();
        FillData(dataEvent,CalibrationRoutines.PassModule);

    }

    public void timerUpdate() {

    }

    public void resetEventListener() {

    }

    public static List<ALERTDetector> getData(DataEvent event){
        ArrayList<ALERTDetector> ALERTList = new ArrayList<ALERTDetector>();
        DataBank adcBank = event.getBank("ATOF::adc");
        DataBank tdcBank = event.getBank("ATOF::tdc");

        ALERTDetector DataSummary = new ALERTDetector();
        DataSummary.setGeometry(
                adcBank.getInt("sector",0),
                adcBank.getInt("superlayer",0),
                adcBank.getInt("layer",0)
        );

        DataSummary.set_TDC(tdcBank.getDouble("tdc_front",0),tdcBank.getDouble("tdc_back",0) );
        DataSummary.setPOS(adcBank.getDouble("xpos",0),
                adcBank.getDouble("ypos",0),
                adcBank.getDouble("zpos",0));
        DataSummary.set_ADC(adcBank.getDouble("adc_front",0),
                adcBank.getDouble("adc_back",0));
        DataSummary.Test_Energy = adcBank.getDouble("energy_front",0);
        DataSummary.YAmp = ALERTDataTransformer(DataSummary).YAmp;
        ALERTList.add(DataSummary);
        return ALERTList;
    }

    public static void FillData(DataEvent event,String name) {
        int sector = 0;
        int super_layer = 0;
        int layer = 0;

        List<ALERTDetector> paddle = new ArrayList<ALERTDetector>();

        ALERTDataStructs Pass = new ALERTDataStructs();
        paddle = getData(event);


        if(name.equals("Pedestal")){
            for (ALERTDetector dr : paddle) {
                super_layer = dr.super_layer;
                sector = dr.sector;
                layer = dr.layer;
                Pedestal[sector][super_layer][layer].fill(dr.ADC_Front);
            }

            if (event.getType() == DataEventType.EVENT_STOP) {
                for (int i =0;i<16;i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k <= 4 - 1; k++) {
                            DataGroup TempGroup2 = new DataGroup(1, 1);
                            TempGroup2.addDataSet(Pedestal[i][j][k], 0);
                            dataGroups.add(TempGroup2, i, j, k);
                        }
                    }
                }
                System.out.println("Passing to DG_Passer");
                Pass.DG_Passer(dataGroups,name);
            }
        }



        else if(name.equals("FBAlign")){
            for (ALERTDetector dr : paddle) {
                super_layer = dr.super_layer;
                sector = dr.sector;
                layer = dr.layer;
                // ** This is being filled with a junk place holder!
                //fbAlign[sector][super_layer][layer].fill(dr.ADC_Front);
                fbAlign[sector][super_layer][layer].fill(dr.HalfTime());
                //Pass.DG_Passer(dataGroups,name);
            }

            if (event.getType() == DataEventType.EVENT_STOP) {
                for (int i =0;i<16;i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k <= 4 - 1; k++) {
                            DataGroup TempGroup2 = new DataGroup(1, 1);
                            TempGroup2.addDataSet(fbAlign[i][j][k], 0);
                            dataGroups.add(TempGroup2, i, j, k);
                        }
                    }
                }
                System.out.println("Passing to DG_Passer");
                Pass.DG_Passer(dataGroups,name);
            }
        }



        else if (name.equals("Veff")) {
            for (ALERTDetector dr : paddle) {
                sector = dr.sector;
                super_layer = dr.super_layer;
                layer = dr.layer;
                //Veff[sector][super_layer][layer].fill(dr.Test_Time,dr.YAmp);
                Veff[sector][super_layer][layer].fill(dr.HalfTime(),dr.YAmp);
                //Veff[sector][super_layer][layer].fill(dr.YAmp, dr.HalfTime());
                System.out.println(event.getType());

            }
            if (event.getType() == DataEventType.EVENT_STOP) {
                System.out.println("the last event");
                for (int i = 0; i < 16; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k <= 4 - 1; k++) {
                            DataGroup TempGroup2 = new DataGroup(1, 1);
                            TempGroup2.addDataSet(Veff[i][j][k], 0);
                            dataGroups.add(TempGroup2, i, j, k);
                        }
                    }
                }
                Pass.DG_Passer(dataGroups, name);
            }
        }

        else if (name.equals("TW")) {
            for (ALERTDetector dr : paddle) {
                super_layer = dr.super_layer;
                sector = dr.sector;
                layer = dr.layer;
                //TW[sector][super_layer][layer].fill(dr.ADC_Front,dr.timeFront);
                TW[sector][super_layer][layer].fill(dr.ADC_Front,dr.HalfTime());
            }

            if (event.getType() == DataEventType.EVENT_STOP) {
                for (int i =0;i<16;i++){
                    for (int j=0;j<2;j++){
                        for (int k=0;k<=4-1;k++) {
                            DataGroup TempGroup2 = new DataGroup(1, 1);
                            TempGroup2.addDataSet(TW[i][j][k], 0);
                            dataGroups.add(TempGroup2, i, j, k);
                        }
                    }
                }
                Pass.DG_Passer(dataGroups,name);
            }
        }
    }

        public void DG_Passer(IndexedList<DataGroup>  ILdg, String name){
            //SCDG = ILdg;
            ALERTCalibrationEngine Send_To_Engine = new ALERTCalibrationEngine();
            Send_To_Engine.CalibHandler(ILdg,name);

            //return SCDG;
        }

    public static ALERTDetector ALERTDataTransformer(ALERTDetector DataSum){
        // apply method in ALERTDetector to get the final form of the data (apply corrections)
        double test_Veff = 0.0;

        ALERTDetector Transformer = new ALERTDetector();
        Transformer.YAmp= Transformer.VeffTestMethod(DataSum.XPOS,DataSum.YPOS,DataSum.Test_ADC);
        //Transformer.VeffTestMethod(ls.XPOS, ls.YPOS));
        return Transformer;
    }





    public void Create_Fill_Histo2D(String Module){
        if (Module.equals("Pedestal")) {
            System.out.println("Initialize Pedestal Histograms");
            for (int i=0; i<16;i++){
                for (int j = 0; j < 2; j++){
                    for (int k=0;k<=4-1;k++){
                        String Hist_Name = String.format("Pedestal_%d_%d_%d", i, j,k);
                        Pedestal[i][j][k] = new H1F("Pedestal",Hist_Name,1000,0,10000);
                    }
                }
            }
        }


       else if (Module.equals("FBAlign")) {
            System.out.println("Initialize FBAlign Histograms");
            for (int i=0; i<16;i++){
                for (int j = 0; j < 2; j++){
                    for (int k=0;k<=4-1;k++){
                        String Hist_Name = String.format("FrontBack_%d_%d_%d", i, j,k);
                        fbAlign[i][j][k] = new H1F("fbAlignment",Hist_Name,20,-10,10);
                    }
                }
            }
        }
        else if (Module.equals("TW")) {
            System.out.println("Initialize Timewalk Histograms");
            for (int i = 0; i < 16; i++) {
                for (int j = 0; j <= 2; j++) {
                    for (int k=0;k<=4-1;k++) {
                        String Hist_Name = String.format("TW_%d_%d_%d", i, j,k);
                        TW[i][j][k]= new H2F("TW", Hist_Name, 500, 10, 2200, 10, -2, 2);
                    }
                }
            }
        }
        else if (Module.equals("Veff")) {
            System.out.println("Initialize Veff Histograms");
            for (int i = 0; i < 16; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k=0;k<=4-1;k++) {
                        String Hist_Name = String.format("Veff_%d_%d_%d", i, j,k);
                        //System.out.println(Hist_Name);
                        //Veff[i][j][k] = new H2F("Veff", Hist_Name, 80, -30, 60, 80, -20, 50);
                        Veff[i][j][k] = new H2F("Veff", Hist_Name, 20, -10, 10, 40, -30, 30);
                    }
                }
            }
        }
    }



}
