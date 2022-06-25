/*
 * This program presents the calculation of:
 *  - The muon self-energy
 *  - The muon anomalous magnetic moment
 *
 *  These calculations are done in a simple QED model 
 *  containing only the muon and the photon, built from
 *  scratch at the beginning of the program.
 *
 *  The two calculations involve loop diagrams and we 
 *  voluntarily show more details than necessary to 
 *  make the reader a bit more used to the MARTY
 *  framework.
 *
 *  Step by step we derive and comment the results, 
 *  how to interpret and use them. 
 *
 *  A numerical library is generated at the end to test
 *  the values calculated at the loop-level by MARTY. 
 *  The script making these numerical tests is 
 *  example_demolib.cpp and must be placed in demolib/script
 *  after the program generation, compiled and executed
 *  (the executable will be in demolib/bin/ after the 
 *  library compilation).
 *
 *  Along the program we refer to specific sections of 
 *  the manual for the methods used:
 *  https://marty.in2p3.fr/doc/marty-manual.pdf
 */
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
    
    // For the general recipee for model building in MARTY see
    // section 5.1

    // Initialize the gauge group
    // See section 5.2 of the manual
    Model model;
    model.addGaugedGroup(group::Type::U1, "em", constant_s("e"));
    model.init();

    model.renameParticle("A_em", "A"); // see section 5.5

    // Create the muon particle
    // See sections 2.1 and 2.2 for the creation of particles
    // and their settings (representation, mass etc).
    Particle muon = diracfermion_s("mu ; \\mu", model);
    muon->setGroupRep("em", -1); // Charge -1 electromagnetic
    muon->setMass(constant_s("m_mu")); 
    model.addParticle(muon);

    // Refresh the model
    model.refresh();

    // Look at what you've done :)
    Display(model); // Model in the terminal

    // For the calculation and interpretation of Feynman rules
    // see section 6.2
    Show(model.getFeynmanRules()); // Feynman diagrams for the vertices

    cout << "Press enter to launch the calculation of the"
              << " muon self-energy ...\n";
    cin.get();

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Calculation of the muon self-energy
    /////////////////////////////////////////////
    /////////////////////////////////////////////

    cout << "###############################\n";
    cout << "####  MUON SELF-ENERGY\n";
    cout << "###############################\n\n";

    // First calculate the one-loop amplitude
    // We take off-shell muon to prevent the
    // Dirac equation to be applied!
    // More details on amplitude calculations in section 6.4
    // Section 6.4.1 explains in particular how to define
    // the particle insertions.
    Amplitude selfEnergy = model.computeAmplitude(
            OneLoop,
            {Incoming(OffShell("mu")), Outgoing(OffShell("mu"))}
            );
    cout << "AMPLITUDE RESULTS:\n";
    Display(selfEnergy); // Amplitude in the terminal
    Show(selfEnergy);    // Feynman diagrams

    // Decompose the amplitude over Lorentz structures
    // The amplitude is multiplied by i automatically
    // Section 6.7 presents in details the extraction of
    // Wilson coefficients in MARTY
    cout << "WILSON COEFFICIENT RESULTS:\n";
    WilsonSet wilsonsSelfEnergy = model.getWilsonCoefficients(selfEnergy); 
    Display(wilsonsSelfEnergy);

    // For each term in the WilsonSet one can obtain the expression of the
    // coefficient using .coef.getCoefficient()
    // (See the "General Wilson coefficient extraction 1/2" sample code
    // in section 6.7.2)
    // The operator expression could be obtained in a similar way using
    // .op.getExpression()
    //
    // Th self-energy contains to terms, one proportional to m_mu
    // and one proportional to \shashed{p} with p the muon momentum.
    // We obtain the two corresponding coefficients in he following.
    Expr muonSelfEnergy_mTerm = wilsonsSelfEnergy[0].coef.getCoefficient();
    Expr muonSelfEnergy_pTerm = wilsonsSelfEnergy[1].coef.getCoefficient();

    // We evaluate the abbreviations to see the exact expression. The list
    // of abbreviations used by MARTY can alo be displayed at any time 
    // using 
    // DisplayAbbreviations();
    cout << "DECOMPOSITION OF THE TWO CONTRIBUTIONS:\n";
    cout << "M-term contribution: " 
              << Evaluated(muonSelfEnergy_mTerm, eval::abbreviation) 
              << endl;
    cout << "P-term contribution: " 
              << Evaluated(muonSelfEnergy_pTerm, eval::abbreviation) 
              << endl;
    cout << endl;

    // We can also compute the squared amplitude if we want
    // See section 6.5 for the calculation of squared amplitudes
    Expr squaredSelfEnergy = model.computeSquaredAmplitude(selfEnergy); 
    cout << "SQUARED AMPLITUDE RESULT:\n";
    // Evaluate the abbreviations
    Expr evaluatedSelfEnergy = Evaluated(squaredSelfEnergy, eval::abbreviation);
    // Simplify by expanding and factoring again
    // As explained below, this is not recommended in general (for large expressions
    // in particular)
    Expr simplifiedSelfEnergy = DeepHardFactored(DeepExpanded(evaluatedSelfEnergy));
    cout << "\nM2              = " << squaredSelfEnergy << endl;
    cout << "\nM2 [evaluated]  = " << evaluatedSelfEnergy << endl;
    cout << "\nM2 [simplified] = " << simplifiedSelfEnergy << endl;

    cout << "\nPress enter to launch the calculation of (g-2) ...\n";
    cin.get();

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Calculation of the muon anomalous
    //  magnetic moment (g-2)
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    
    cout << "###############################\n";
    cout << "####  MUON MAGNETIC MOMENT\n";
    cout << "###############################\n\n";

    // Here we can directly compute the Wilson coefficients 
    // as we do not square the amplitude
    WilsonSet wilsonsMuonVertex = model.computeWilsonCoefficients(
            OneLoop,
            {Incoming("mu"), Outgoing("mu"), Outgoing("A")}
            );
    cout << "WILSON COEFFICIENTS RESULTS:\n";
    Display(wilsonsMuonVertex);
    Show(wilsonsMuonVertex);

    // To get the contribution of a particular operator, we first
    // need to create the operator. 
    // In this case the (chromo-)magnetic operator for the muon
    // (mu*sigma_{mu,nu}*mu)*F^{mu,nu}
    // The definition of operators for the extraction of Wilson
    // coefficients is exhaustively presented in section 6.7.2
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

    cout << "MAGNETIC MOMENT RESULTS:";
    cout << "Muon magnetic moment              = "
              << Evaluated(muonMagneticMoment, eval::abbreviation)
              << endl;
    // Remove abbreviations using Evaluated(expr, eval::abbreviation)
    Expr evaluatedMagneticMoment = Evaluated(muonMagneticMoment, eval::abbreviation);
    cout << "Muon magnetic moment [evaluated]  = "
              << evaluatedMagneticMoment
              << endl;

    // Simplify small expressions with DeepHardFactored(DeepExpanded())
    // This is however not recommended on large expressions!
    // For pedagocical purposes and on small results this is however really good :)
    Expr simplifiedMagneticMoment = DeepHardFactored(DeepExpanded(evaluatedMagneticMoment));
    cout << "Muon magnetic moment [simplified] = "
              << simplifiedMagneticMoment
              << endl;

    cout << "\nPress enter to launch the library generation ...\n";
    cin.get();

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //  Generation of the library
    //  (this is the easy part :D)
    /////////////////////////////////////////////
    /////////////////////////////////////////////

    // For the code generation and to use the generated library
    // see the chapter 7 of the manual :))
    Library lib("demolib");

    // In case we change function names it is better to remove old files 
    // in the library to ensure the new library does not contain old and
    // inconsistent files:
    lib.cleanExistingSources();  

    // We add the functions one by one, giving only the name and symbolic
    // expression to compile
    lib.addFunction("mu_self_e_mterm", muonSelfEnergy_mTerm);
    lib.addFunction("mu_self_e_pterm", muonSelfEnergy_pTerm);
    lib.addFunction("mu_self_e_squared", squaredSelfEnergy);
    lib.addFunction("mu_magnetic_vertex", muonMagneticMoment);
    lib.addFunction("mu_magnetic_vertex_eval", evaluatedMagneticMoment);
    lib.addFunction("mu_magnetic_vertex_simpli", simplifiedMagneticMoment);

    // Make MARTY build automatically the library :)
    // We could also use a simple
    // lib.print();
    // if we want to compile th library later on (useful for large libraries)
    lib.build();


    return 0;
}
