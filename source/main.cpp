#include <iostream>
#include <regex>
#include "mfr.hpp"

int main(int argc, char *argv[]) {
    std::vector<std::string> args(argv, argv + argc);
    std::string inputFile;
    if (argc == 2) {
       inputFile = args[1];
    }else{
        std::cout << "You have to add the command line argument for the input file." << std::endl;
        return 0;
    }
    MicroFlowReactor::WeakFlame mfr(inputFile);

    auto sol = newSolution(mfr.MechanismFile, "gas");
    auto gas = sol->thermo();

    MicroFlowReactor::WeakFlameStFlow flow(gas, mfr.NusseltNumber, mfr.TubeDiameter);

    const auto NumSpecies = gas->nSpecies();
    const auto SpeciesNames = gas->speciesNames();

    //TODO: We have to discuss how to replace ",".
    std::ofstream outfile(mfr.FileNameOut, std::ios::trunc);
    outfile << "x (m),Temperature (K),Tw(K),Velocity (m/s),Heat Release Rate (J/s/m3)", Density (kg/m3), Pressure (MPa)";
    for (auto i = 0; i < NumSpecies; i++) {
        std::string replacedSpecName = std::regex_replace(SpeciesNames[i], std::regex(","),"__");
        outfile << ",X_" << replacedSpecName;
    }
    outfile << "\n";

    vector_fp x(NumSpecies, 0.0);
    if (mfr.SimWithPhi) {
        gas->setEquivalenceRatio(mfr.EquivalenceRatio, mfr.Fuel, mfr.Oxidizer);
        gas->setState_TP(mfr.Temperature, mfr.Pressure);
    } else {
        gas->setState_TPX(mfr.Temperature, mfr.Pressure, mfr.Mixture);
    }

    gas->getMoleFractions(x.data());
    mfr.Rho = gas->density();
    mfr.MassFlowRate = mfr.Rho * mfr.Velocity;

    flow.readTwFile(mfr.FileNameIn);

    flow.setAxisymmetricFlow();

    vector_fp grid = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0};
    const auto NumGrid = grid.size();
    vector_fp z(NumGrid), T(NumGrid), u(NumGrid);
    for (int iz = 0; iz < NumGrid; iz++) {
        z[iz] = grid[iz] * mfr.TubeLength;
        T[iz] = flow.getTwAt(z[iz]);
        u[iz] = T[iz] / mfr.Temperature * mfr.Velocity;
    }
    flow.setupGrid(NumGrid, &z[0]);

    std::unique_ptr<Transport> tran(newTransportMgr("Mix", sol->thermo().get()));
    flow.setTransport(*tran);
    flow.setKinetics(*sol->kinetics());
    flow.setPressure(mfr.Pressure);
    flow.setSteadyTolerances(mfr.RtolMRF, mfr.AtolMRF);

    Inlet1D inlet;
    inlet.setMoleFractions(x.data());
    inlet.setMdot(mfr.MassFlowRate);
    inlet.setTemperature(mfr.Temperature);

    Outlet1D outlet;

    std::vector<Domain1D *> domains{&inlet, &flow, &outlet};
    Sim1D flame(domains);

    if (mfr.isRestart){
        flame.restore(mfr.FileNameRestart, "solution");
    } else{
        flame.setInitialGuess("velocity", grid, u);
        flame.setInitialGuess("T", grid, T);
    }
    flow.solveEnergyEqn();
    flame.showSolution();
    flame.setRefineCriteria(mfr.Domain, mfr.RefineCriteriaRatio, mfr.RefineCriteriaCurve, mfr.RefineCriteriaSlope);
    inlet.setMoleFractions(x.data());
    inlet.setMdot(mfr.MassFlowRate);
    inlet.setTemperature(mfr.Temperature);

    flame.solve(mfr.Log, true);
    flame.writeStats();

    std::cout << "===== Converged =====" << std::endl;

    // Postprocessing
    // TODO: Create some functions in WeakFlame class for postprocessing.
    for (auto n = 0; n < flow.nPoints(); n++) {
        // gas conditions
        auto tmpTemperature = flame.value(mfr.Domain, flow.componentIndex("T"), n);
        auto tmpPressure = mfr.Pressure;
        vector_fp tmpMassFrac(NumSpecies, 0.0);
        vector_fp tmpMoleFrac(NumSpecies, 0.0);
        for (auto i = 0; i < NumSpecies; i++) tmpMassFrac[i] = flame.value(1, flow.componentIndex(SpeciesNames[i]), n);

        gas->setState_TP(tmpTemperature, tmpPressure);
        gas->setMassFractions(&tmpMassFrac[0]);

        auto kin = sol->kinetics();
        vector_fp wdot(NumSpecies, 0.0);
        kin->getNetProductionRates(&wdot[0]);
        vector_fp hk(NumSpecies, 0.0);
        gas->getPartialMolarEnthalpies(&hk[0]);
        auto HeatReleaseRate = 0.0;
        for (auto i = 0; i < NumSpecies; i++) HeatReleaseRate += -wdot[i] * hk[i];
        gas->getMoleFractions(tmpMoleFrac.data());

        outfile << flow.grid(n) << ","
                << tmpTemperature << ","
                << flow.getTwAt(flow.grid(n)) << ","
                << flame.value(mfr.Domain, flow.componentIndex("velocity"), n) << ","
                << HeatReleaseRate << ","
                << flow.density(n) << ","
                << flow.pressure()/1e6; // MPa
        for (auto i = 0; i < NumSpecies; i++) outfile << "," << tmpMoleFrac[i];
        outfile << "\n";
    }

    if (mfr.WriteRawData) {
        std::cout << "Save xml file.";
        flame.save(mfr.FileNameRaw, "solution", "");
    }
}
