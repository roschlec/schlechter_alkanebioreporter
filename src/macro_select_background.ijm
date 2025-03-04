// Set global measurements
run("Fresh Start");
run("Set Measurements...", "area mean modal min display redirect=None decimal=3");

// Set the directory of the images and an output directory
#@ File (label="Select input directory", style = "directory") dir
#@ File (label="Select output directory", style = "directory") out

list = getFileList(dir); // Get a list of all files to analyze

setBatchMode(true);

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".czi")) {
        print("Analysing image " + list[i]);
        openfile = "open=[" + dir + "/" + list[i] + "] autoscale=false";
        run("Bio-Formats Windowless Importer", openfile);

        // Convert to 8-bit grayscale if necessary
        run("8-bit");

        // Get names and remove extension
        titleID = getTitle();
        atitle = replace(titleID, ".czi", "");

        // Duplicate for segmentation
        run("Duplicate...", "duplicate channels=1 title=Mask");

        // Apply threshold based on title to detect cells
        setAutoThreshold("Percentile dark");

        run("Convert to Mask");
        run("Invert"); // Invert to get background

        // Create a selection of the background area
        run("Create Selection");

        // Switch back to the original image
        selectWindow(titleID);

        // Set ROI Manager
        roiManager("reset");

        // Number of random background squares
        numROIs = 50;
        squareSize = 2;

        // Get image dimensions
        imageWidth = getWidth();
        imageHeight = getHeight();

        count = 0;
        while (count < numROIs) {
            x = random() * (imageWidth - squareSize);
            y = random() * (imageHeight - squareSize);

            makeRectangle(x, y, squareSize, squareSize);

            // Ensure ROI is in the background (outside cells)
            if (selectionType() != -1) {
                roiManager("Add");
                count++;
            }
        }

        // Measure background intensity
        roiManager("Measure");

        // Save results
        saveAs("Results", out + "/" + atitle + "_background.csv");

        // Cleanup
        roiManager("reset");
        run("Clear Results");
        close(); // Close original image
        selectWindow("Mask");
        close(); // Close mask image
    }
}

setBatchMode(false); // End batch mode