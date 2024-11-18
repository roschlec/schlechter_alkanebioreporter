// Convert to 8-bit grayscale if necessary
run("8-bit");

// Get the histogram of the image
getStatistics(area, mean, min, max, std);
thr = round(mean - 4*std);
setThreshold(0, thr);

// Analyze particles
run("Set Measurements...", "area mean perimeter fit shape feret's display redirect=None decimal=3");
run("Analyze Particles...", "size=0.40-2.50 exclude display add");
resetThreshold();

// Get the image title
imageTitle = getTitle();  // This retrieves the current image's title

// Get the number of ROIs in the ROI Manager
n = roiManager("count");

// Create a new Results table (if not already open)
if (isOpen("Results")) {
    run("Clear Results");
} else {
    run("Results...");
}

roiManager("Show None");
// Loop through each detected ROI (particle)
for (i = 0; i < n; i++) {
    roiManager("select", i);  // Select the current ROI

    run("Measure");

    // Ask the user to classify this particle as a bacterial cell using getBoolean()
    message = "Is particle " + (i + 1) + " a bacterial cell?";  // Custom message for each particle
    if (getBoolean(message, "Yes", "No")) {
        cellStatus = "Yes";  // User confirmed it is a bacterial cell
    } else {
        cellStatus = "No";   // User denied it
    }

    // Add the result to the table (ROI number and cell status)
    setResult("ROI", i, i + 1);  // ROI number (1-based indexing)
    setResult("Cell Status", i, cellStatus);  // User's classification (Yes/No)

    updateResults();  // Update the Results table

}

// Remove the ROI after measurement
roiManager("deselect");
roiManager("delete");

// Save the results with the image title (ensure the file name doesn't contain invalid characters)
safeFileName = replace(imageTitle, " ", "_");  // Replace spaces with underscores
safeFileName = replace(safeFileName, ".", "_"); // Replace periods with underscores (avoid file extension issues)

// Save the results as a CSV file
saveAs("Results", "/Users/rschlechter/repo/schlechter_alkanebioreporter/data/train_data/" + safeFileName + "_train.csv");

close();