//	Process alkane bioreporter images
	
// 	Set global measurements
run("Fresh Start");
run("Set Measurements...", "area mean perimeter fit shape feret's display redirect=None decimal=3");

//	Set the directory of the images and an output directory
#@ File (label="Select input directory", style = "directory") dir
#@ File (label= "Select ouput directory", style = "directory") out

list = getFileList(dir); // Get a list of all the files that will be analysed

setBatchMode(true);
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], ".czi")){
		print("Analysing image " + list[i]);
		openfile = "open=[" + dir + "/" + list[i] + "] autoscale=false";
		run("Bio-Formats Windowless Importer", openfile);
		//open(dir+"/"+list[i]);

		// Convert to 8-bit grayscale if necessary
		run("8-bit");
		
		// Get names and split 
		titleID = getTitle();
		atitle = replace(titleID, ".czi", "");
		run("Split Channels");

		// Phase contrast
		selectWindow("C2-" + titleID);
		pc = getTitle();
		pcID = getImageID();

		selectWindow("C1-" + titleID);
		gfp = getTitle();
		gfpID = getImageID();
		
		// THRESHOLD the PHASE CONTRAST image
		selectImage(pcID);
		getStatistics(area, mean, min, max, std);
		if (indexOf(titleID, "lb") >= 0) {
			thr = round(mean - 1*std);
			} else if (indexOf(titleID, "diesel") >= 0) {
				thr = round(mean - 2*std);
				} else {
					thr = round(mean - 4*std);
					}
		setThreshold(0, thr);

		// Create a mask
		selectImage(pcID);
		run("Convert to Mask");
		run("Erode");
		run("Erode");
		run("Dilate");
		run("Dilate");
		run("Watershed");

		//	MEASURE GFP IMAGE
		selectImage(pcID);
		run("Analyze Particles...", "size=0.20-2.50 exclude display add");
		
		// Check if any ROIs were detected
    	roiManager("Count");
    	roiCount = roiManager("Count");

    	// Create a new Results table (if not already open)
		if (isOpen("Results")) {
		    run("Clear Results");
		} else {
		    run("Results...");
		}
    
    	if (roiCount > 0) {
        	// Save the ROIs to a file
        	selectImage(gfpID);
			roiManager("Measure");

			// Remove the ROI after measurement
			roiManager("deselect");
			roiManager("delete");
    	} else {
        	print("No ROIs detected for: " + atitle + ", skipping...");
        }
        
        // Save results table for the image
    	saveAs("Results", out + "/" + atitle + "_results.csv");
        
		}
		while (nImages()>0) {
          selectImage(nImages());  
          run("Close");
    }
}