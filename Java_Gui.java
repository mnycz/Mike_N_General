package org.clas.gui;

import org.clas.modules.ALERTCalibrationEngine;
import org.clas.modules.ALERTDataStructs;
import org.clas.modules.geom.ALERTGeometry;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.CalibrationConstantsListener;
import org.jlab.detector.calib.utils.CalibrationConstantsView;
import org.jlab.detector.view.DetectorListener;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedList;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class ALERTCalibGUI implements IDataEventListener, ActionListener, CalibrationConstantsListener, DetectorListener{

    private final int npaddles = 4;
    //public IndexedList<DataGroup> dataGroups_Constants = new IndexedList<DataGroup>(3);
    //public ALERTCalConstants ce = null;
    public JPanel mainPanel = null;
    public DataSourceProcessorPane processorPane = null;
    public JSplitPane splitPanel = null;
    public JPanel detectorPanel = null;
    public DetectorPane2D detectorView = null;
    public JSplitPane moduleView = null;
    public EmbeddedCanvas canvas = null;
    public CalibrationConstantsView ccview = null;
    public ALERTCalibrationEngine ce = null;


    public ALERTCalibGUI(){
        ce = new ALERTCalibrationEngine();
        ccview = new CalibrationConstantsView();

        ccview.addConstants(ce.getCalibrationConstants().get(0));

        mainPanel = new JPanel();
        mainPanel.setLayout(new BorderLayout());

        // create detector panel
        detectorPanel = new JPanel();
        detectorPanel.setLayout(new BorderLayout());
        detectorView = new DetectorPane2D();
        drawDetector();
        detectorView.getView().addDetectorListener(this);
        detectorPanel.add(detectorView);


        moduleView = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        canvas = new EmbeddedCanvas();
        ccview = new CalibrationConstantsView();
        //System.out.println("GetCalCon"+ce.getCalibrationConstants());


        ccview.addConstants(ce.getCalibrationConstants().get(0),this);
        moduleView.setTopComponent(canvas);
        moduleView.setBottomComponent(ccview);
        moduleView.setDividerLocation(0.5);
        moduleView.setResizeWeight(0.6);

        splitPanel = new JSplitPane();
        //splitPanel.setLeftComponent(null);
        splitPanel.setLeftComponent(detectorView);
        splitPanel.setRightComponent(moduleView);
        processorPane = new DataSourceProcessorPane();
        processorPane.setUpdateRate(10000);
        processorPane.addEventListener(this);



        mainPanel.add(splitPanel);
        mainPanel.add(processorPane, BorderLayout.PAGE_END);
    }

    public void actionPerformed(ActionEvent e) {

    }

    public void constantsEvent(CalibrationConstants calibrationConstants, int i, int i1) {

        System.out.println("WHEN DOES CONSTANTS Events get called");
        //ALERTCalibrationEngine ce2 = new ALERTCalibrationEngine();
        System.out.println("Well. it's working " + i + "  " + i1);
        String str_sector    = (String) calibrationConstants.getValueAt(i1, 0);
        String str_layer     = (String) calibrationConstants.getValueAt(i1, 1);
        String str_component = (String) calibrationConstants.getValueAt(i1, 2);
        System.out.println(str_sector + " " + str_layer + " " + str_component);
        IndexedList<DataGroup> group = ce.getDataGroup();

        int sector    = Integer.parseInt(str_sector);
        int layer     = Integer.parseInt(str_layer);
        int component = Integer.parseInt(str_component);

        
        DataGroup dataGroup = group.getItem(sector,layer,component);
       
        this.canvas.clear();
        this.canvas.draw(dataGroup);
        this.canvas.update();
        resetEventListener();


    }


    public void dataEventAction(DataEvent dataEvent) {
        //System.out.println("In Calib GUI");
        ALERTCalibrationEngine CalibrationRoutines = new ALERTCalibrationEngine();
        ALERTDataStructs Passing = new ALERTDataStructs();
        Passing.FillData(dataEvent,  CalibrationRoutines.PassModule);

    }

    public void timerUpdate() {

    }

    public void resetEventListener() {

    }


    public void drawDetector() {
        //System.out.println("Draw Detector called");
        double FTOFSize = 500.0;
        int[] widths = new int[] {4,4,4};
        int[] lengths =new int[] {2,2,2};

        //String[]  names    = new String[]{"FTOF 1A","FTOF 1B","FTOF 2"};
        String[] names = new String[]{"ALERT_ATOF_1","ALERT_ATOF_2","ALERT_ATOF_3"};
        //for(int sector = 1; sector <= 5; sector++){
        for(int sector = 1; sector <= 15; sector++){
            //double rotation = Math.toRadians((sector-1)*(360.0/6)+90.0);
            double rotation = Math.toRadians((sector-1)*(360.0/15)+90.0);
            for(int layer = 1; layer <=2; layer++){
                int width  = widths[layer-1];
                int length = lengths[layer-1];
                for(int paddle = 1; paddle <= 4; paddle++){
                    ALERTGeometry BAR_XZ = new ALERTGeometry();
                    DetectorShape2D shape = new DetectorShape2D();
                    //shape.getDescriptor().setType(DetectorType.FTOF);
                    shape.getDescriptor().setType(DetectorType.UNDEFINED);
                    shape.getDescriptor().setSectorLayerComponent(sector, layer, paddle);
                    //shape.createBarXY(10 + length*paddle, width);
                    shape.createBarXY(1 + length*2.0, width);
                    shape.getShapePath().translateXYZ(0.0, 40 + width*paddle , 0.0);
                    //shape.getShapePath().translateXYZ(0.0, 40 + width*paddle , 0);
                    shape.getShapePath().rotateZ(rotation);
                 
                    detectorView.getView().addShape(names[layer-1], shape);
                   
                }
            }
        }
        detectorView.setName("ATOF");
        detectorView.updateBox();
    }




    public void processShape(DetectorShape2D dsd) {
        // show summary
        int sector = dsd.getDescriptor().getSector();
        int layer  =  dsd.getDescriptor().getLayer();
        int paddle = dsd.getDescriptor().getComponent();
        System.out.println("Selected shape " + sector + " " + layer + " " + paddle);
        IndexedList<DataGroup> group = ce.getDataGroup();

        if(group.hasItem(sector,layer,paddle)==true){
            this.canvas.clear();
            this.canvas.draw(this.ce.getDataGroup().getItem(sector,layer,paddle));
            this.canvas.update();
        } else {
            System.out.println(" ERROR: can not find the data group");
        }

    }




    public void Draw(H2F blank, F1D fit) {
// This is drawing an addititonal panel for each Calibration Method. Should be tweaked in order
// to use the same Frame? An additional button to move through the different calibration methods?
        System.out.println("IN main");
        JFrame frame = new JFrame("Calibration");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ALERTCalibGUI viewer = new ALERTCalibGUI();
        viewer.canvas.draw(blank);
        viewer.canvas.draw(fit,"same");

        frame.add(viewer.mainPanel);
        frame.setSize(1400, 800);
        frame.setVisible(true);
    }

    public void Draw_H1(H1F blank, F1D fit) {
// This is drawing an addititonal panel for each Calibration Method. Should be tweaked in order
// to use the same Frame? An additional button to move through the different calibration methods?
        System.out.println("IN main");
        JFrame frame = new JFrame("Calibration");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ALERTCalibGUI viewer = new ALERTCalibGUI();
        viewer.canvas.draw(blank);
        viewer.canvas.draw(fit, "same");

        frame.add(viewer.mainPanel);
        frame.setSize(1400, 800);
        frame.setVisible(true);

    }

}
