//
// Created by Youhi Morii on 2021/06/21.
//

#ifndef MFR_MFR_HPP
#define MFR_MFR_HPP

#include <fstream>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/kinetics.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/Solution.h"
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/transport.h"

using namespace Cantera;

namespace MicroFlowReactor {

    struct WeakFlame {
        std::string Fuel;
        std::string Oxidizer;
        std::string Mixture;
        double EquivalenceRatio;
        double Velocity;
        double Pressure;
        double Temperature;

        double TubeLength;
        double TubeDiameter;

        // Weak flame
        double NusseltNumber;
        std::string FileNameIn;
        std::string FileNameOut;
        std::string FileNameRaw;
        std::string FileNameRestart;
        std::string MechanismFile;

        const int Domain = 1;
        int Log;

        double RefineCriteriaRatio;
        double RefineCriteriaSlope;
        double RefineCriteriaCurve;

        double Rho;
        double MassFlowRate;
        bool WriteRawData;
        bool SimWithPhi;
        bool isRestart;

        double AtolMRF;
        double RtolMRF;

        WeakFlame(std::string InputFile) {
            try {
                YAML::Node config = YAML::LoadFile(InputFile);
                // PrePost block
                WriteRawData = config["PrePost"]["rawData"].as<bool>();
                FileNameRaw = config["PrePost"]["rawFile"].as<std::string>("rawData.yaml");
                // cantera block
                MechanismFile = config["cantera"]["mechFile"].as<std::string>();
                // mfr block
                SimWithPhi = config["mfr"]["simWithPhi"].as<bool>();
                isRestart = config["mfr"]["Restart"].as<bool>();
                FileNameRestart = config["mfr"]["restartFile"].as<std::string>();
                if (SimWithPhi) {
                    EquivalenceRatio = config["mfr"]["equilRatio"].as<double>();
                    Fuel = config["mfr"]["fuelComp"].as<std::string>();
                    Oxidizer = config["mfr"]["oxyComp"].as<std::string>();
                } else {
                    Mixture = config["mfr"]["mixComp"].as<std::string>();
                }
                Velocity = config["mfr"]["uVel"].as<double>();
                Pressure = config["mfr"]["pres"].as<double>() * Cantera::OneAtm;
                Temperature = config["mfr"]["temp"].as<double>();
                TubeLength = config["mfr"]["length"].as<double>();
                TubeDiameter = config["mfr"]["diameter"].as<double>();
                FileNameIn = config["mfr"]["tempWall"].as<std::string>();
                FileNameOut = config["mfr"]["outFile"].as<std::string>();
                NusseltNumber = config["mfr"]["nusseltNum"].as<double>();
                AtolMRF = config["mfr"]["atol"].as<double>();
                RtolMRF = config["mfr"]["rtol"].as<double>();
                // solver block
                Log = config["solver"]["logging"].as<int>();
                // refineCriterion block
                RefineCriteriaRatio = config["refineCriterion"]["ratio"].as<double>();
                RefineCriteriaSlope = config["refineCriterion"]["slope"].as<double>();
                RefineCriteriaCurve = config["refineCriterion"]["curve"].as<double>();
            }
            catch (YAML::Exception &e) {
                std::cerr << e.what() << std::endl;
                exit(-1);
            }
        };

        ~WeakFlame() {};
    };

    struct WeakFlameStFlow : public Cantera::StFlow {
        double NusseltNumber;
        double TubeDiameter; //[m]

        WeakFlameStFlow(auto gas, double nusseltNumber, double tubeDiameter) : Cantera::StFlow(gas) {
            NusseltNumber = nusseltNumber;
            TubeDiameter = tubeDiameter;
        };

        ~WeakFlameStFlow() {};

        void readTwFile(std::string fileName);

        double getTwAt(double positionX);

        void evalResidual(double *x, double *rsd, int *diag, double rdt, size_t jmin, size_t jmax) override;

    private:
        std::vector<double> PositionZ;
        std::vector<double> TemperatureWall;

        // TODO: add ROP
        // TODO: Hdot
    };

    // TODO: Add class for FREI simulation
    struct FREI {
    };

}

#endif //MFR_MFR_HPP
