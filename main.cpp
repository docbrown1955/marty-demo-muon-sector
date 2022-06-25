#include "marty.h"

using namespace std;
using namespace csl;
using namespace mty;

int main() 
{
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Model definition
    /////////////////////////////////////////////
    /////////////////////////////////////////////

    // Initialize the gauge group
    Model model;
    model.addGaugedGroup(group::Type::U1, "em", constant_s("e"));
    model.init();

    model.renameParticle("A_em", "A");

    // Create the muon particle
    Particle muon = diracfermion_s("mu ; \\mu", model);
    muon->setGroupRep("em", -1); // Charge -1 electromagnetic
    muon->setMass(constant_s("m_mu")); 
    model.addParticle(muon);

    // Refresh the model
    model.refresh();

    // Look at what you've done :)
    Display(model); // Model in the terminal
    Show(model.getFeynmanRules()); // Feynman diagrams for the vertices

    std::cout << "Press enter to launch the calculation of the"
              << " muon self-energy ...\n";
    std::cin.get();

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Calculation of the muon self-energy
    /////////////////////////////////////////////
    /////////////////////////////////////////////

    std::cout << "###############################\n";
    std::cout << "####  MUON SELF-ENERGY\n";
    std::cout << "###############################\n\n";

    // First calculate the one-loop amplitude
    // We take off-shell muon to prevent the
    // Dirac equation to be applied!
    Amplitude selfEnergy = model.computeAmplitude(
            OneLoop,
            {Incoming(OffShell("mu")), Outgoing(OffShell("mu"))}
            );
    std::cout << "AMPLITUDE RESULTS:\n";
    Display(selfEnergy); // Amplitude in the terminal
    Show(selfEnergy);    // Feynman diagrams

    // Decompose the amplitude over Lorentz structures
    // The amplitude is multiplied by i
    std::cout << "WILSON COEFFICIENT RESULTS:\n";
    WilsonSet wilsonsSelfEnergy = model.getWilsonCoefficients(selfEnergy); 
    Display(wilsonsSelfEnergy);

    Expr muonSelfEnergy_mTerm = wilsonsSelfEnergy[0].coef.getCoefficient();
    Expr muonSelfEnergy_pTerm = wilsonsSelfEnergy[1].coef.getCoefficient();

    std::cout << "DECOMPOSITION OF THE TWO CONTRIBUTIONS:\n";
    std::cout << "M-term contribution: " 
              << Evaluated(muonSelfEnergy_mTerm, eval::abbreviation) 
              << std::endl;
    std::cout << "P-term contribution: " 
              << Evaluated(muonSelfEnergy_pTerm, eval::abbreviation) 
              << std::endl;
    std::cout << std::endl;

    // We can also compute the squared amplitude if we want
    Expr squaredSelfEnergy = model.computeSquaredAmplitude(selfEnergy); 
    // To print an expression without abbreviation use the following syntax:
    std::cout << "SQUARED AMPLITUDE RESULT:\n";
    Expr evaluatedSelfEnergy = Evaluated(squaredSelfEnergy, eval::abbreviation);
    Expr simplifiedSelfEnergy = DeepHardFactored(DeepExpanded(evaluatedSelfEnergy));
    std::cout << "\nM2              = " << squaredSelfEnergy << std::endl;
    std::cout << "\nM2 [evaluated]  = " << evaluatedSelfEnergy << std::endl;
    std::cout << "\nM2 [simplified] = " << simplifiedSelfEnergy << std::endl;

    std::cout << "\nPress enter to launch the calculation of (g-2) ...\n";
    std::cin.get();

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Calculation of the muon anomalous
    //  magnetic moment (g-2)
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    
    std::cout << "###############################\n";
    std::cout << "####  MUON MAGNETIC MOMENT\n";
    std::cout << "###############################\n\n";

    // Here we can directly compute the Wilson coefficients 
    // as we do not square the amplitude
    WilsonSet wilsonsMuonVertex = model.computeWilsonCoefficients(
            OneLoop,
            {Incoming("mu"), Outgoing("mu"), Outgoing("A")}
            );
    std::cout << "WILSON COEFFICIENTS RESULTS:\n";
    Display(wilsonsMuonVertex);
    Show(wilsonsMuonVertex);

    // To get the contribution of a particular operator, we first
    // need to create the operator. 
    // In this case the (chromo-)magnetic operator for the muon
    // (mu*sigma_{mu,nu}*mu)*F^{mu,nu}
    vector<Wilson> muonMagOp = chromoMagneticOperator(
            model, 
            wilsonsMuonVertex, 
            DiracCoupling::S  // Scalar coupling after sigma_mu_nu
    );
    // Finally we extract the coefficient of the particular operator
    // we received from chromoMagneticOperator()
    Expr muonMagneticMoment = getWilsonCoefficient(
            wilsonsMuonVertex, 
            muonMagOp
    );

    std::cout << "MAGNETIC MOMENT RESULTS:";
    std::cout << "Muon magnetic moment              = "
              << Evaluated(muonMagneticMoment, eval::abbreviation)
              << std::endl;
    // Remove abbreviations using Evaluated(expr, eval::abbreviation)
    Expr evaluatedMagneticMoment = Evaluated(muonMagneticMoment, eval::abbreviation);
    std::cout << "Muon magnetic moment [evaluated]  = "
              << evaluatedMagneticMoment
              << std::endl;

    // Simplify small expressions with DeepHardFactored(DeepExpanded())
    // This is however not recommended on large expressions!
    // For pedagocical purposes on small results this is however really good :)
    Expr simplifiedMagneticMoment = DeepHardFactored(DeepExpanded(evaluatedMagneticMoment));
    std::cout << "Muon magnetic moment [simplified] = "
              << simplifiedMagneticMoment
              << std::endl;

    std::cout << "\nPress enter to launch the library generation ...\n";
    std::cin.get();

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Generation of the library
    //  (this is the easy part :D)
    /////////////////////////////////////////////
    /////////////////////////////////////////////

    Library lib("demolib");
    lib.cleanExistingSources();
    lib.addFunction("mu_self_e_mterm", muonSelfEnergy_mTerm);
    lib.addFunction("mu_self_e_pterm", muonSelfEnergy_pTerm);
    lib.addFunction("mu_self_e_squared", squaredSelfEnergy);
    lib.addFunction("mu_magnetic_vertex", muonMagneticMoment);
    lib.addFunction("mu_magnetic_vertex_eval", evaluatedMagneticMoment);
    lib.addFunction("mu_magnetic_vertex_simpli", simplifiedMagneticMoment);
    lib.build();


    return 0;
}
