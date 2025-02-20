//	Process alkane bioreporter images - get background values
	
// 	Set global measurements
run("Fresh Start");
run("Set Measurements...", "mean display redirect=None decimal=3");

//	Set the directory of the images and an output directory
#@ File (label="Select input directory", style = "directory") dir
#@ File (label= "Select ouput directory", style = "directory") out

list = getFileList(dir); // Get a list of all the files that will be analysed

setBatchMode(true);

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], '.czi')) {
        print("Analysing image " + list[i]);
        openfile = "open=[" + dir + "/" + list[i] + "]";
        run("Bio-Formats Windowless Importer", openfile);

        // Convert to 8-bit grayscale if necessary
        run("8-bit");

        // Get image dimensions
        imageWidth = getWidth();
        imageHeight = getHeight();

        // Get names and split
        titleID = getTitle();
        atitle = replace(titleID, ".czi", "");

        // Number of random points to generate
        n = 8;

        // Generate random points across the image
        for (j = 0; j < n; j++) {
            
            // Generate random coordinates within the image dimensions
            x = random() * imageWidth;
            y = random() * imageHeight;

            // Create a random point ROI at the generated coordinates
            makePoint(x, y);
            roiManager("Add"); // Add the point ROI to the ROI Manager
        }

        // Measure the intensity of the points in the ROI Manager
        roiManager("measure");
        roiManager("deselect");
        roiManager("delete");

        // Save results table for the image
    	saveAs("Results", out + "/" + atitle + "_background.csv");
        Table.deleteRows(0, n-1);

		}
		while (nImages()>0) {
          selectImage(nImages());  
          run("Close");
    }
}

setBatchMode(false); // End batch mode