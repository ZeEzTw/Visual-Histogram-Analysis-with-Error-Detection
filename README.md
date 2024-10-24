# Visual Histogram Analysis with Error Detection

This project provides tools for analyzing histograms from ROOT files and detecting potential issues based on predefined criteria. The main goal is to visually highlight problematic areas and generate a JSON report for further analysis.

## Features

- **Histogram Visualization**: Renders histograms from ROOT files with customizable overlays.
- **Error Detection**: Identifies issues by checking if data points fall outside specified regions.
- **Visual Overlays**: Adds circles and boxes to highlight problematic areas on histograms.
- **JSON Reporting**: Generates a JSON file with details of any detected issues for further analysis.

## How to Run It

1. Open ROOT with the command:
   
        root --web=off
 
2. Load and execute the code using:

       .L histogram_checker.C

## File Descriptions

- `histogram_checker.C`: Main script for visual histogram analysis and error detection.
- `problems.json`: JSON file containing details of detected issues.
- `output/`: Directory for any output files generated by the script.

##Example JSON Output

The JSON report lists the histograms with issues and specifies which region (circle or rectangle) was violated. We have two rectangles: a red rectangle and circle (upper bound) and a green rectangle and circle (lower bound).

#Example JSON:

```
{
  "Histogram_with_problems": "hEE_0_2",
  "Issues": [
    "The point is outside the green rectangle",
    "The point is not within the green circle"
  ]
}
```
