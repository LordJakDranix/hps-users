package org.hps.users.mrsolt;

import java.util.Arrays;

import org.lcsim.event.EventHeader;

/**
 * This is the tuple template driver
 * Use this to add your code and variables to make a tuple
 * Run the GeneralTupleDriver to output info into a text file
 * Change the steering file to include this driver
 * Run "makeTree.py" on text file to create a root tuple
 *
 * @author mrsolt on Aug 31, 2017
 */

public class TupleTemplateDriver extends GeneralTupleDriver {
	
    public static void fillVariables(EventHeader event) {
    	
    	//Fill tuple with run number and event number
        tupleMap.put("run/I", (double) event.getRunNumber());
        tupleMap.put("event/I", (double) event.getEventNumber());
        
        
        /* Put code here
         * <some code>*/
        

        /* Fill variables here. For example, to add number of tracks "nTrk", do:
         *  tupleMap.put("nTrk/I", (double) ntrk); 
         *  where "I" means it is an integer*/
    }
    
    public static void addVariables() {
    	/* Add variables to the list "newVars"
    	 * Do it in the form "varName/<varType>" 
    	 * including the quotes and separated by commas
    	 * <varType> is "I" for integer, "D" for double, and "B" for boolean*/
    	
        String[] newVars = new String[]{"run/I", "event/I"};
        tupleVariables.addAll(Arrays.asList(newVars));
    }
}