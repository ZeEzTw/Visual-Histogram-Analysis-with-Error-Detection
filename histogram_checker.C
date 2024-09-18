{
    const int ndet = 32;
    const int nE = 3;
    const double Eg[nE] = {778.9045, 964.072, 1408.013};

    TCanvas *canvas = new TCanvas("c", "Histograms", 800, 600);
    TFile *outFile = new TFile("histograms_with_circles_and_squares.root", "RECREATE");

    gSystem->Exec("mkdir -p output");
    std::ofstream jsonFile("problems.json");

    jsonFile << "[\n";

    std::map<std::pair<int, int>, double> errorMap;
    std::ifstream infile("CL32_matrix.dat");
    int i, j;
    double errorValue;

    // Load error values into the map
    while (infile >> i >> j >> errorValue)
    {
        errorMap[std::make_pair(i, j)] = errorValue;
    }

    bool firstProblem = true;

    // Loop through all detector combinations
    for (int i = 0; i < ndet; i++)
    {
        for (int j = i + 1; j < ndet; j++)
        {
            TFile *fEE = new TFile(Form("hEE_%i_%i.root", i, j));
            if(!fEE || fEE->IsZombie())
            {
                delete fEE;
                continue;
            }
            TH2D *hEE;
            fEE->GetObject(Form("hEE%i %i", i, j), hEE);
            if (!hEE)
            {
                delete fEE;
                continue;
            }

            hEE->GetXaxis()->SetRangeUser(0, 1500);
            hEE->GetYaxis()->SetRangeUser(0, 1500);
            hEE->Draw("COLZ");

            // Draw histograms for each energy
            /*for (int k = 0; k < nE; k++)
            {
                TFile *f = TFile::Open(Form("hEE_%i_%i_%i.root", i, j, k));
                if (!f || f->IsZombie() || !f)
                {
                    delete f;
                    continue; // Skip this iteration
                }
                //TFile *f = new TFile(Form("hEE_%i_%i_%i.root", i, j, k));
                TH2D *h;
                f->GetObject(Form("h_%i_%i_%i", i, j, k), h);
                if (h)
                {
                    h->Draw("same COLZ");
                    delete h;
                }
                f->Close();
                delete f;
            }
            */

            auto key1 = std::make_pair(i + 1, j + 1);
            auto key2 = std::make_pair(j + 1, i + 1);

            if (errorMap.find(key1) != errorMap.end() && errorMap.find(key2) != errorMap.end())
            {
                // Calculate coordinates and draw circles
                double errorValue1 = errorMap[key1];
                double xCoord1 = 1173 + errorValue1 * 1332;
                double yCoord1 = 1332 + errorValue1 * 1173;
                TEllipse *circle1 = new TEllipse(xCoord1, yCoord1, 15, 15);
                circle1->SetLineColor(kRed);
                circle1->SetFillStyle(0);
                circle1->SetLineWidth(2);
                circle1->Draw("same");

                double errorValue2 = errorMap[key2];
                double xCoord2 = 1332 + errorValue2 * 1173;
                double yCoord2 = 1173 + errorValue2 * 1332;
                TEllipse *circle2 = new TEllipse(xCoord2, yCoord2, 15, 15);
                circle2->SetLineColor(kGreen);
                circle2->SetFillStyle(0);
                circle2->SetLineWidth(2);
                circle2->Draw("same");

                double smallSize = 20, largeSize = 40;

                TBox *boxLarge1 = new TBox(xCoord1 - largeSize / 2.0, yCoord1 - largeSize / 2.0,
                                           xCoord1 + largeSize / 2.0, yCoord1 + largeSize / 2.0);
                boxLarge1->SetLineColor(kBlack);
                boxLarge1->SetLineWidth(2);
                boxLarge1->SetFillStyle(0);
                boxLarge1->Draw("same");

                TBox *boxLarge2 = new TBox(xCoord2 - largeSize / 2.0, yCoord2 - largeSize / 2.0,
                                           xCoord2 + largeSize / 2.0, yCoord2 + largeSize / 2.0);
                boxLarge2->SetLineColor(kBlack);
                boxLarge2->SetLineWidth(2);
                boxLarge2->SetFillStyle(0);
                boxLarge2->Draw("same");

                double smallSum1 = 0, largeSum1 = 0, smallSum2 = 0, largeSum2 = 0;

                auto calculateSum = [&](double xMin, double xMax, double yMin, double yMax) -> double
                {
                    double sum = 0;
                    int binXMin = hEE->GetXaxis()->FindBin(xMin);
                    int binXMax = hEE->GetXaxis()->FindBin(xMax);
                    int binYMin = hEE->GetYaxis()->FindBin(yMin);
                    int binYMax = hEE->GetYaxis()->FindBin(yMax);

                    for (int binX = binXMin; binX <= binXMax; ++binX)
                    {
                        for (int binY = binYMin; binY <= binYMax; ++binY)
                        {
                            sum += hEE->GetBinContent(binX, binY);
                        }
                    }
                    return sum;
                };

                smallSum1 = calculateSum(xCoord1 - smallSize / 2.0, xCoord1 + smallSize / 2.0,
                                         yCoord1 - smallSize / 2.0, yCoord1 + smallSize / 2.0);
                largeSum1 = calculateSum(xCoord1 - largeSize / 2.0, xCoord1 + largeSize / 2.0,
                                         yCoord1 - largeSize / 2.0, yCoord1 + largeSize / 2.0);
                smallSum2 = calculateSum(xCoord2 - smallSize / 2.0, xCoord2 + smallSize / 2.0,
                                         yCoord2 - smallSize / 2.0, yCoord2 + smallSize / 2.0);
                largeSum2 = calculateSum(xCoord2 - largeSize / 2.0, xCoord2 + largeSize / 2.0,
                                         yCoord2 - largeSize / 2.0, yCoord2 + largeSize / 2.0);

                largeSum1 -= smallSum1;
                largeSum2 -= smallSum2;

                if (largeSum1 > smallSum1 || largeSum2 > smallSum2 || smallSum1 == 0 || smallSum2 == 0)
                {
                    if (!firstProblem)
                    {
                        jsonFile << ",\n";
                    }
                    firstProblem = false;

                    jsonFile << "{\n";
                    jsonFile << "\t\"Histogram_with_problems\": \"" << Form("hEE_%i_%i", i, j) << "\",\n";

                    bool hasRedCircleIssue = largeSum1 > smallSum1 || smallSum1 == 0;
                    bool hasGreenCircleIssue = largeSum2 > smallSum2 || smallSum2 == 0;

                    jsonFile << "\t\"Issues\": [\n";

                    if (hasRedCircleIssue && hasGreenCircleIssue)
                    {
                        jsonFile << "\t\t\"The point is not in the circle (Red circle) and the point is not in the circle (Green circle)\"\n";
                    }
                    else if (hasRedCircleIssue)
                    {
                        jsonFile << "\t\t\"The point is not in the circle (Red circle)\"\n";
                    }
                    else if (hasGreenCircleIssue)
                    {
                        jsonFile << "\t\t\"The point is not in the circle (Green circle)\"\n";
                    }

                    jsonFile << "\t]\n";
                    jsonFile << "}";
                }

                canvas->Update();
                canvas->SaveAs(Form("output/histogram_%i_%i.root", i, j));
                outFile->WriteTObject(hEE, Form("histogram_%i_%i", i, j));
            }
            delete hEE;
            delete fEE;
        }
    }

    jsonFile << "\n]\n";
    jsonFile.close();

    outFile->Close();
}
