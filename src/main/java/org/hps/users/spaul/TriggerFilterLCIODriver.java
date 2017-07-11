package org.hps.users.spaul;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.hps.recon.filtering.EventReconFilter;
import org.hps.recon.filtering.Single0TriggerFilterDriver;
import org.hps.record.epics.EpicsData;
import org.hps.record.scalers.ScalerData;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TIData;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.lcio.LCIOWriter;
import org.lcsim.util.Driver;
import org.lcsim.util.loop.LCIODriver;
/**
 * hacked version of LCIODriver, which writes out events with the specified trigger.  
 * @author spaul
 *
 */
public class TriggerFilterLCIODriver extends Driver {
    
    
    private LCIOWriter writer;
    private Set<String> listIgnore = new HashSet<String>();
    private Set<String> listKeep = new HashSet<String>();
    private File outputFile;

    private boolean writeSingle0 = false;

    private boolean writeSingle1  = false;

    private boolean writePair1  = false;

    private boolean writePulser  = false;

    private boolean writePair0  = false;
    
    public void setWriteSingle1(boolean b){
        writeSingle1 = b;
    }
    
    public void setWriteSingle0(boolean b){
        writeSingle0 = b;
    }
    
    public void setWritePair0(boolean b){
        writePair0 = b;
    }
    
    public void setWritePair1(boolean b){
        writePair1 = b;
    }
    
    public void setWritePulser(boolean b){
        writePulser = b;
    }

    public TriggerFilterLCIODriver(String file) {
        this(addFileExtension(file), null);
    }

    public TriggerFilterLCIODriver(File file) {
        this(file, null);
    }

    public TriggerFilterLCIODriver(String file, Collection<String> listIgnore) {
        this(new File(addFileExtension(file)), listIgnore);
    }

    public TriggerFilterLCIODriver(File file, Collection<String> listIgnore) {
        this.outputFile = file;
        if (listIgnore != null) {
            this.listIgnore.addAll(listIgnore);
        }
    }

    public TriggerFilterLCIODriver() {
    }

    public void setOutputFilePath(String filePath) {
        outputFile = new File(addFileExtension(filePath));
    }

    public void setIgnoreCollections(String[] ignoreCollections) {
        listIgnore.addAll(Arrays.asList(ignoreCollections));
    }

    public void setWriteOnlyCollections(String[] keepCollections) {
        listKeep.addAll(Arrays.asList(keepCollections));
    }

    public void setIgnoreCollection(String ignoreCollection) {
        listIgnore.add(ignoreCollection);
    }

    public void setWriteOnlyCollection(String writeOnlyCollection) {
        listKeep.add(writeOnlyCollection);
    }

    private void setupWriter() {
        // Cleanup existing writer.
        if (writer != null) {
            try {
                writer.flush();
                writer.close();
                writer = null;
            } catch (IOException x) {
                System.err.println(x.getMessage());
            }
        }

        // Setup new writer.
        try {
            writer = new LCIOWriter(outputFile);
        } catch (IOException x) {
            throw new RuntimeException("Error creating writer", x);
        }
        writer.addAllIgnore(listIgnore);
        writer.addAllWriteOnly(listKeep);

        try {
            writer.reOpen();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    protected void startOfData() {
        setupWriter();
    }

    protected void endOfData() {
        try {
            writer.close();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    protected void process(EventHeader event) {
        
     // 1. keep all events with EPICS data (could also use event tag = 31):
        if (EpicsData.read(event) != null) return;

        // 2. keep all events with Scaler data:
        if (ScalerData.read(event) != null) return;

        // 3. drop event if it doesn't have a TriggerBank
        if (!event.hasCollection(GenericObject.class,"TriggerBank"))
          throw new Driver.NextEventException();
      
        // 4. keep event if it was from the chosen trigger:
        boolean writeThisEvent = false;
        for (GenericObject gob : event.get(GenericObject.class,"TriggerBank"))
        {
          if (!(AbstractIntData.getTag(gob) == TIData.BANK_TAG)) continue;
          TIData tid = new TIData(gob);
          if (tid.isSingle0Trigger() && writeSingle0 || 
                  tid.isSingle1Trigger() && writeSingle1 || 
                  tid.isPair1Trigger() && writePair1 || 
                  tid.isPulserTrigger() && writePulser ||
                  tid.isPair0Trigger() && writePair0) writeThisEvent = true;
        }
        if(! writeThisEvent)
            return;
        
        try {
            writer.write(event);
        } catch (IOException x) {
            throw new RuntimeException("Error writing LCIO file", x);
        }
    }

    protected void suspend() {
        try {
            writer.flush();
        } catch (IOException x) {
            throw new RuntimeException("Error flushing LCIO file", x);
        }
    }

    private static String addFileExtension(String filePath) {
        if (!filePath.endsWith(".slcio")) {
            return filePath + ".slcio";
        } else
            return filePath;
    }
    
}
