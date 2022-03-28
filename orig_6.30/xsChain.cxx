#include <XSstreams.h>
#include <xsTypes.h>
#include <XSContainer.h>
#include <XSFit/Fit/Fit.h>
#include <XSFit/MCMC/ChainManager.h>
#include <XSFit/MCMC/AsciiChain.h>
#include <XSFit/MCMC/FITSChain.h>
#include <XSFit/Randomizer/RandomizerBase.h>
#include <XSModel/Model/Model.h>
#include <XSUser/Global/Global.h>
#include <XSUser/Global/XSGlobal.h>
#include <XSUser/Handler/HandlerUtils.h>
#include <XSUser/Handler/XSinterface.h>
#include <XSUtil/Parse/XSparse.h>
#include <XSUtil/Signals/SignalHandler.h>
#include <XSUtil/Utils/XSutility.h>
#include <memory>
#include <sstream>
#include <set>

int
XSGlobal::xsChain(ClientData cdata,Tcl_Interp* tclInterp,int objc, Tcl_Obj* CONST  objv[])
{
   StringArray rawArgs;
   HandlerUtils::tclArgsToCpp(objc, objv, rawArgs);
   int status = doChain(rawArgs);
   if (!status)
      return globalData->autoSave(TCL_OK);
   else
      return globalData->autoSave(TCL_ERROR);
}

int XSGlobal::doChain(const StringArray & rawArgs)
{
  enum SUBCOMS {BEST=1, BURN, CLEAR, DIC, FILETYPE, INFO, LENGTH, LOAD, 
		PROPOSAL, RAND, RECALC, RESCALE, RUN, STAT, TEMPERATURE, TYPE, 
		UNLOAD, WALKERS};
  static std::map<string, size_t> subComs;
  static size_t burn = 0;
  static size_t length = 0;
  static size_t walkers = 10;
  static Real temperature = 1.0;
  static bool randomize = false;
  static string parSpec("1");
  static RangePair prevRange(0,0);
  static string fileType(FITSChain::FORMATNAME());
  static string chainType("gw");
  const size_t nArgs = rawArgs.size();
  string cmd("chain");
  string subcom;
  using namespace XSContainer;

  if (!subComs.size()) {
    subComs["best"] = (size_t)BEST;
    subComs["burn"] = (size_t)BURN;
    subComs["clear"] = (size_t)CLEAR;
    subComs["dic"] = (size_t)DIC;
    subComs["filetype"] = (size_t)FILETYPE;
    subComs["info"] = (size_t)INFO;
    subComs["length"] = (size_t)LENGTH;
    subComs["load"] = (size_t)LOAD;
    subComs["proposal"] = (size_t)PROPOSAL;
    subComs["rand"] = (size_t)RAND;
    subComs["recalc"] = (size_t)RECALC; // Aliased to proposal
    subComs["rescale"] = (size_t)RESCALE;
    subComs["run"] = (size_t)RUN;
    subComs["stat"] = (size_t)STAT;
    subComs["temperature"] = (size_t)TEMPERATURE;
    subComs["type"] = (size_t)TYPE;
    subComs["unload"] = (size_t)UNLOAD;
    subComs["walkers"] = (size_t)WALKERS;
  }
  if (nArgs == 1 || (subcom = rawArgs[1]) == "?") {
    XSutility::printValidOptions(tcout, cmd, subComs);
  } else {
    subcom = XSutility::lowerCase(subcom);
    std::map<string,size_t>::const_iterator itSubCom(subComs.lower_bound(subcom));
    std::map<string,size_t>::const_iterator itSubComEnd(subComs.end());
    size_t choice = 0;
    if (itSubCom != itSubComEnd && itSubCom->first.find(subcom) == 0)
      choice = itSubCom->second;
    string thirdArg = (nArgs > 2) ? rawArgs[2] : string("");
    string fourthArg = (nArgs > 3) ? rawArgs[3] : string("");
    std::istringstream ssThirdArg(thirdArg);
    ChainManager* chainManager = fit->chainManager();
    const ChainManager::ChainContainer& loadedChains = chainManager->chains();
    const size_t width = chainManager->width();
    try {
      size_t tmpInput=0;
      string tmpMsg;
      size_t nChains = loadedChains.size();
      
      switch (choice) {
      case BEST:
	// write out the statistic value and best-fit parameters in the loaded
	// chains
	if (!nChains) {
	  tcout << "\nNo loaded chains" << endl;
	} else {
	  RealArray parVals;
	  Real statVal;
	  chainManager->findBestFit(parVals, statVal);
          // All chains have the same parameters.  Let's just get the
          //  labels from the first.
          const Chain* firstChain = loadedChains.begin()->second;
          const vector<Chain::ParamID>& parIDs = firstChain->paramIDs();
          
	  tcout << "\nLowest fit statistic value found: " << statVal << endl;
	  tcout << "    for parameter values : " << endl;
          for (size_t i=0; i<parVals.size(); ++i)
          {
             const Chain::ParamID& parID = parIDs[i];
             tcout.width(15);
             if (parID.modName != Model::DEFAULT())
                tcout << parID.modName << ":";
             tcout << parID.index << "   ";
	     tcout.width(10);
	     tcout << parVals[i] << endl;             
          }
          
	}
	break;
      case BURN:
      case LENGTH:
	tmpMsg = (choice == BURN) ? string("burn"):string("length");
	if (nArgs > 2) {
	  if (!(ssThirdArg >> tmpInput) || !ssThirdArg.eof()) {
	    tcerr << "   Integer length parameter required: chain "
		  << tmpMsg << " <length>" << std::endl;
	    return -1;
	  } else {
	    choice == BURN ? burn = tmpInput : length = tmpInput;
	  }
	} else {
	  tcout << "   Integer length parameter required: chain "
		<< tmpMsg << " <length>" << std::endl;
	}
	break;
      case CLEAR:
	// chainManager will automatically set its
	// isSynchedWithFitParams flag to false.
	chainManager->clearChains();
	tcout << "   All chains are now removed." << std::endl;
	break;
      case DIC:
	// write out the Deviance Information Criterion and effective number
	// of parameters from the loaded chains
	Real effNumPars, devInfCrit;
	chainManager->calcDevInfCrit(effNumPars, devInfCrit);
	tcout << "Deviance Information Criterion: " << devInfCrit << endl;
	tcout << "Effective number of Parameters: " << effNumPars << endl;
	break;
      case FILETYPE:
	if (nArgs > 2) {
	  string fType = rawArgs[2];
	  fType = XSutility::lowerCase(fType);
	  if (AsciiChain::FORMATNAME().find(fType) == 0) {
	    tcout << "New chain files will be in ASCII text format."
		  << std::endl;
	    fileType = AsciiChain::FORMATNAME();
	  } else if (FITSChain::FORMATNAME().find(fType) == 0) {
	    tcout << "New chain files will be in FITS format."
		  << std::endl;
	    fileType = FITSChain::FORMATNAME();
	    
	    //  tcout << "FITS chain format is not implemented yet." << std::endl;
	  } else if ( fType == "?" ) {
	    tcout << "       Valid formats are: "
		  <<  AsciiChain::FORMATNAME() << "|" << FITSChain::FORMATNAME()
		  << "\n       Current setting is:  " << fileType << std::endl;
	    return -1;      
	  } else {
	    tcerr << "***Error:  Unrecognized chain file format type."
		  << "\n       Valid formats are: "
		  <<  AsciiChain::FORMATNAME() << "|" << FITSChain::FORMATNAME()
		  << "\n       Current setting is:  " << fileType << std::endl;
	    return -1;      
	  }
	} else {
	  tcout << "Choose output format: chain filetype "
		<< AsciiChain::FORMATNAME() << "|" << FITSChain::FORMATNAME()
		<<"\n      Current setting is:  " << fileType << std::endl;
	}
	break;
      case INFO:
	{
	  using namespace std;
	  tmpMsg = randomize ? string("on"):string("off");
	  tcout << "Output file format           : " << fileType
		<< "\nChain type                   : " << chainType 
		<< "\nCurrent chain length setting : " << length 
		<< "\nBurn-in length               : " << burn
		<< "\nChain width                  : " << width;
	  if ( chainType == "gw" ) tcout << "\nNumber of walkers            : " << walkers;
	  if ( chainType == "mh" ) {
	    tcout << "\nCurrent temperature          : ";
	    int savePrec = tcout.precision(Chain::TEMPERPREC());
	    tcout << showpoint << scientific << temperature << noshowpoint;
	    tcout.setf(ios_base::fmtflags(0), ios_base::floatfield);
	    tcout.precision(savePrec);
	    tcout << "\nRandomization " << tmpMsg << endl;
	    chainManager->reportChainProposal();
	  }
	  if (!nChains) {
	    tcout << "\nNo loaded chains" << endl;
	  } else {
	    tcout << "\nLoaded chains:   length:" << endl;
	    ChainManager::ChainContainer::const_iterator it =
	      loadedChains.begin();
	    ChainManager::ChainContainer::const_iterator itEnd = 
	      loadedChains.end();
	    size_t index=1;
	    while (it != itEnd) {
	      tcout << "    " << index << ".  " << it->first 
		    << "  " << it->second->length() << endl;
	      ++index;
	      ++it;
	    }
	  }
	}
	break;
      case LOAD:
	if (nArgs < 3) {
	  tcout << "   File name parameter required: chain load <fileName>"
		<< std::endl;
	} else {
	  if (chainManager->getChain(thirdArg)) {
	    string prompt(thirdArg);
	    prompt += " already exists as a loaded chain.";
	    prompt += "\nDo you want to replace it? (y/n): ";
	    string result;
	    XSparse::basicPrompt(prompt, result);
	    if (result.empty() || (result[0] != 'y'&&
				   result[0] != 'Y'))  break;
	  }
	  // Both the ChainIO and Chain contructors may throw when
	  // creating from an existing file.
	  std::unique_ptr<ChainIO> apIO((ChainIO*)0);
	  // We have no idea what format this file is in or even if
	  // it exists, so...
	  try {
	    apIO.reset(new AsciiChain(thirdArg));
	  } catch (ChainIO::ChainIOError&) {
	    // If in here, file cannot even be opened.
	    throw;
	  } catch (YellowAlert&) {
	    // File exists, but not in the proper text format.
	    try {
	      apIO.reset(new FITSChain(thirdArg));
	    } catch (YellowAlert&) {
	      tcerr << "***Error: Chain file has an unrecognized format."
		    << std::endl;
	      throw;
	    }
	  }
	  std::unique_ptr<Chain> apChain(new Chain(apIO.get(),thirdArg));
	  // If Chain ctor succeeded, it will now own ChainIO*.
	  apIO.release();
	  // addChain is guaranteed to add chain to
	  // container or else throw. 
	  chainManager->addChain(apChain.get());
	  // chainManager now owns new chain.
	  apChain.release();
	  if (chainManager->chains().size() == 1) {
	    length = chainManager->length();
	    tcout << "  Default chain length is now set to length of "
		  << thirdArg << ": " << length << std::endl;
	    const Chain* firstChain=chainManager->chains().begin()->second;
	    if (firstChain->chainType()=="GoodmanWeare") {
	      if (walkers != firstChain->walkers()) {
		walkers = firstChain->walkers();
		tcout << "  Default walkers length is now set to: "
		      << walkers << std::endl;
	      }
	    }                                                 
	    chainManager->isSynchedWithFitParams(fit->checkChainsForSynch());
	  }
	}
	break;
      case RECALC: // Now aliased to "chain proposal <dist> chain"
      case PROPOSAL:
	if (nArgs < 3 && choice != RECALC) {
	  tcout << "   Currently loaded chain proposal options: "
		<< fit->getChainProposalNames() << std::endl;
	} else {
	  string propName;
	  size_t optionalArgsStart = string::npos;
	  size_t optionalArgsEnd = nArgs;
	  bool isNative = true;
	  string initString;
	  const string WS(" \t\n");
	  if (choice == RECALC) {
	    const string& currentName = 
	      chainManager->chainProposal()->name();
	    string::size_type spaceLoc = currentName.find_first_of(WS);
	    // Only built-in proposals can have whitespace
	    if (spaceLoc != string::npos) {
	      string distName = currentName.substr(0,spaceLoc);
	      propName = distName + " chain";
	    } else {
	      tcerr << "  Current chain proposal does not use covariance matrix."
		    << "\n  No chain proposal recalculation will be performed."
		    << std::endl;
	      return -1;                        
	    }
	  } else {
	    propName = XSutility::lowerCase(thirdArg);
	    if (propName[0] == '?') {
	      tcout << "   Currently loaded chain proposal options: "
		    << fit->getChainProposalNames() << std::endl;
	      return 0;                        
	    }
	    // propName may be an abbreviation at this point, but
	    // getRandomizingStrategy can deal with it.

	    const RandomizerBase* testProp = fit->getRandomizingStrategy(propName);
	    if (!testProp) {
	      string errMsg("Unrecognized chain proposal name: ");
	      errMsg += propName;
	      errMsg += "\nValid proposal names: ";
	      errMsg += XSContainer::fit->getChainProposalNames();
	      throw Chain::ChainError(errMsg);
	    }

	    string::size_type wsLoc = testProp->name().find_first_of(WS);
	    propName = testProp->name().substr(0,wsLoc);
	    // propName now is the FULL name of the first word
	    // of testProp's name.

	    const std::set<string>& nativeNames = fit->nativeRandomizerNames();
	    isNative = nativeNames.find(propName) != nativeNames.end();
	    if (isNative) {
	      // Currently existing built-in classes make no use of optional
	      // strings except for the case where the 4th arg is to be 
	      // interpreted as a filename.
	      if (nArgs == 3) {
		tcerr << "***Error: To use " << propName << " proposal, must also specify "
		      << "\nfit | chain | diagonal | matrix | deltas | <filename>" << std::endl;
		return -1;
	      }
	      const string sourceQualifier(XSutility::lowerCase(fourthArg));

	      // sourceQualifier may be an abbreviation for fit|chain|diagonal|matrix|deltas.
	      // (We're going to assume this is more likely than someone trying to 
	      // name a file "fi" or "cha".)
	      if (string("fit").find(sourceQualifier) == 0) propName += " fit";
	      else if (string("chain").find(sourceQualifier) == 0) propName += " chain";
	      else if (string("diagonal").find(sourceQualifier) == 0) {
		propName += " <cmdline>";
		initString = "diagonal ";
		optionalArgsStart = 4;
	      } else if (string("matrix").find(sourceQualifier) == 0) {
		propName += " <cmdline>";
		initString = "matrix ";
		optionalArgsStart = 4;
	      } else if (string("deltas").find(sourceQualifier) == 0) {
		propName += " deltas";
		initString = " ";
		optionalArgsStart = 4;
	      } else {
		// Assume fourthArg is a filename, to be passed in initString.
		optionalArgsStart = 3;
		optionalArgsEnd = 4;
		// This label is needed to match file-handling 
		// Randomizer template classes.
		propName += " <filename>";
	      }
	    } else { // end if native
	      // For user add-on classes, the name can only be
	      // specified with the 3 arg (hence it must be 
	      // 1 word, no spaces).
	      optionalArgsStart = 3;
	    }
	  }  // end if !="chain recalc"
	  for (size_t i=optionalArgsStart; i<optionalArgsEnd; ++i) {
	    initString += rawArgs[i];
	    if (i < optionalArgsEnd-1) initString += " ";
	  }
	  
	  chainManager->setChainProposal(propName, initString);
	  if (isNative && propName.find(" chain") != string::npos) chainManager->recalc();
                     
	  if (chainType == "gw") {
	    tcout << "    Note that the currently selected chain type is Goodman-Weare so the the 'chain proposal'"
		  <<" setting is only used to generate the initial walkers." << std::endl;
	  }
	} // end if nArgs >= 3
	break;
      case RAND:
	thirdArg = XSutility::lowerCase(thirdArg);
	if (thirdArg == "on") {
	  randomize = true;
	} else if (thirdArg == "off") {
	  randomize = false;
	} else {
	  tcout << "   Must specify chain rand <on|off>" << std::endl;
	}
	if (chainType == "gw") {
	  tcout << "    Note that rand has no effect on Goodman-Weare." << std::endl;
	}
	break;
      case RESCALE:
	{
	  Real factor = 0.0;
	  if (!(ssThirdArg >> factor) || !ssThirdArg.eof()) {
	    tcout << "   Rescale factor required: chain rescale <factor>"
		  << std::endl;
	  } else {
	    if ( factor == 0.0 ) {
	      tcerr << "  Rescale factor cannot be 0.0." << std::endl;
	      return -1;
	    }
	    // try to get the current covariance matrix
	    const RealArray* pCovarMatrix = 
	      fit->chainManager()->chainProposal()->covariance(fit);
	    if ( pCovarMatrix ) {
	      const RealArray& covarMatrix = *pCovarMatrix;
	      StringArray propParts;
	      string propName = fit->chainManager()->chainProposal()->name();
	      XSparse::collateByWhitespace(propParts, propName);
	      stringstream ssInit;
	      ssInit << "matrix ";
	      // we need the lower half and diagonal of the covariance matrix
	      const size_t nPar = static_cast<size_t>
		(sqrt(static_cast<Real>(covarMatrix.size())));
	      for (size_t i=0; i<nPar; i++) {
		for (size_t j=0; j<=i; j++) {
		  ssInit << factor * covarMatrix[i*nPar+j] << " ";
		}
	      }
	      string initString = ssInit.str();
	      chainManager->setChainProposal(propParts[0]+" <cmdline>", initString);
	    }
	  }
	}
	break;
      case RUN:
	if (nArgs < 3) {
	  tcout << "   Must specify a file name: chain run <fileName>" 
		<< std::endl;
	} else {
	  if (!length) {
	    string msg("Chain length is zero - please set it.");
	    throw Chain::ChainError(msg);
	  }
                  
	  // Need to check for ">" specifier before name,
	  // indicating append, or "!" specifier indicating
          // that it's coming from PyXspec with append=False.
          // In the latter case we want to overwrite any existing
          // file without prompting.
	  string fileName(thirdArg+fourthArg);
	  bool isAppend = false;
          bool isPyXspecOverwrite = false;
	  if (fileName[0] == '>') {
	    isAppend = true;
	    fileName = fileName.substr(1);
	    if (!fileName.length()) {
	      tcout << "  Must specify a file name: chain run <fileName>" 
		    <<std::endl;
	      break;
	    }
	  }
          else if (fileName[0] == '!') {
             isPyXspecOverwrite = true;
             fileName = fileName.substr(1);
	     if (!fileName.length()) {
	       tcout << "  Must specify a file name: chain run <fileName>" 
		     <<std::endl;
	       break;
	     }             
          }
          

	  if (isAppend) {
	    const Chain* chain = chainManager->getChain(fileName);
	    if (!chain) {
	      throw Chain::ChainError("Chain must already be loaded to append to it.");
	    }
	    if (chain->chainType()=="GoodmanWeare") {
	      if (chainType != "gw") {
		string err("To append to a Goodman-Weare chain, chain type must be set to 'gw'");
		throw Chain::ChainError(err);
	      }
	      if (walkers != chain->walkers()) {
		ostringstream oss;
		oss << "This GW chain may only be appended to with walkers = "
		    << chain->walkers();
		throw Chain::ChainError(oss.str());
	      }
	    }
	    if (chain->chainType()=="MetropolisHastings" && chainType != "mh") {
	      string err("To append to a Metropolis-Hastings chain, chain type must be set to 'mh'");
	      throw Chain::ChainError(err);
	    }
	  }
                  
	  // Constraints for gw only: length and burn must be an
	  //  integer multiple of walkers, and walkers must be even.
	  if (chainType == "gw") {
	    if (!walkers) {
	      string msg("Chain number of walkers is zero - please set it.");
	      throw Chain::ChainError(msg);
	    }
	    if (walkers > length) {
	      string msg("Chain length is currently less than number of walkers.\n");
	      msg +="It must be set to a multiple of the number of walkers.";
	      throw Chain::ChainError(msg);
	    }
	    if ( length % walkers != 0 ) {
	      length -= (length % walkers);
	      tcout << "  Length must be a multiple of the number of walkers." << std::endl;
	      tcout << "  Modifying length to " << length << std::endl;
	    }
	    if ( burn % walkers != 0 ) {
	      burn -= (burn % walkers);
	      tcout << "  Burn-in length must be a multiple of the number of walkers." << std::endl;
	      tcout << "  Modifying burn-in to " << burn << std::endl;
	    }                     
	  }
                  
	  if (isAppend) {
	    SIGINT_Handler intHandler;
	    SignalHandler* sigContainer = SignalHandler::instance();
	    EventHandler* oldHandler = sigContainer->registerHandler(SIGINT, &intHandler);
	    try {
	      if ( chainType == "mh" ) {
		chainManager->appendToChain(fileName, length, temperature, fileType);
	      } else if ( chainType == "gw" ) {
		chainManager->appendToChain(fileName, length, walkers, fileType);
	      }
	    } catch (...) {
	      sigContainer->registerHandler(SIGINT,oldHandler);
	      throw;
	    }
	    sigContainer->registerHandler(SIGINT,oldHandler);
	  } else {
	    // Check for pre-existing chain with same file name.
	    // If it exists, must prompt user (unless PyXspec's
            // run() call contained append=False).
	    if (!isPyXspecOverwrite && chainManager->getChain(fileName)) {
	      string prompt(fileName);
	      prompt += " already exists as a loaded chain.";
	      prompt += "\nDo you want to replace it? (y/n): ";
	      string result;
	      XSparse::basicPrompt(prompt, result);
	      if (result.empty() || (result[0] != 'y'&&
				     result[0] != 'Y'))  break;

	    } else if (!isPyXspecOverwrite) {
	      // Chain with filename is not loaded. 
	      // Now check if filename actually exists as a file. 
	      std::ifstream testFile(fileName.c_str());
	      if (testFile) {
		string prompt(fileName);
		prompt += " already exists.  Overwrite? (y/n): ";
		string answer;
		XSparse::basicPrompt(prompt, answer);
		if (!answer.length() || (answer[0] != 'y' && answer[0] != 'Y'))
		  break;
	      }
	    }
	    ChainIO* pIO = 0;
	    // These constructors won't throw when creating a brand
	    // new chain.
	    if (fileType == AsciiChain::FORMATNAME()) pIO = new AsciiChain();
	    else // Assume this must be FITS, there are no other choices.
	      pIO = new FITSChain();
	    // Chain object will now own ChainIO and is wholly responsible
	    // for its destruction.
	    std::unique_ptr<Chain> apChain((Chain*)0);
	    if ( chainType == "mh" ) {
	      apChain.reset(new Chain(pIO, fileName, burn, length, randomize));
	    } else if ( chainType == "gw" ) {
	      apChain.reset(new Chain(pIO, fileName, burn, length, walkers));
	    }
	    SIGINT_Handler intHandler;
	    SignalHandler* sigContainer = SignalHandler::instance();
	    EventHandler* oldHandler = sigContainer->registerHandler(SIGINT,
								     &intHandler);
	    try {
	      apChain->run(0, temperature);
	    } catch (Chain::ChainInterrupt&) {
	      // If chain was interrupted, it will not be loaded.
	      // We also don't want to leave an earlier chain loaded 
	      // that happens to have the same filename.  
	      // The user would have expected that to be overwritten
	      // anyway.
	      sigContainer->registerHandler(SIGINT,oldHandler);
	      chainManager->removeChain(fileName);
	      // The single removal removeChain function does not 
	      // automatically call checkLengths like removeChains
	      // and addChain do (to avoid needless calls when it
	      // is run many consecutive times).  So call it here...
	      chainManager->checkLengths();
	      throw;
	    } catch (...) {
	      sigContainer->registerHandler(SIGINT,oldHandler);
	      throw;
	    }
	    sigContainer->registerHandler(SIGINT,oldHandler);
	    // addChain may throw
	    chainManager->addChain(apChain.get());
	    apChain.release();
	    if (chainManager->chains().size()==1)
	      chainManager->isSynchedWithFitParams(fit->checkChainsForSynch());
	  } // end if not append
	}
	break;
      case STAT:
	if (nArgs > 2) {
	  string modName;
	  size_t iPar=0;
	  if ((convertParameterLabel(thirdArg, modName, iPar)
	       && iPar > 0)|| iPar > 0) {
	    chainManager->calcStat(iPar, modName);
	    parSpec = thirdArg;
	  } else {
	    tcerr << thirdArg << " is not a valid parameter specifier."
		  << std::endl;
	    return -1;
	  }
	} else {
	  // parSpec will always have proper syntax,
	  // no need to check.
	  string modName;
	  size_t iPar=0;
	  convertParameterLabel(parSpec, modName, iPar);
	  chainManager->calcStat(iPar, modName);
	}
	break;
      case TEMPERATURE:
	if (nArgs < 3) {
	  tcout << "  Temperature value required: chain temperature <temp val>"
		<< std::endl;
	} else {
	  string temperStr(rawArgs[2]);
	  std::istringstream iss(temperStr);
	  Real tempTest = 0.0;
	  if (!(iss >> tempTest) || !iss.eof()) {
	    tcerr << "  " << temperStr << " is not a valid temperature value."
		  << std::endl;
	    return -1;
	  }
	  if (tempTest == 0.0) {
	    tcerr << "  Temperature = 0.0 is not allowed." << std::endl;
	    return -1;
	  }
	  temperature = tempTest;
	}
	break;
      case TYPE:
	if (nArgs > 2) {
	  string type = rawArgs[2];
	  type = XSutility::lowerCase(type);
	  if (type == "mh") {
	    tcout << "New chain runs will use Metropolis-Hastings." << std::endl;
	    chainType = "mh";
	  } else if (type == "gw") {
	    tcout << "New chain runs will use Goodman-Weare." << std::endl;
	    chainType = "gw";
	  } else if ( type == "?" ) {
	    tcout << "       Valid types are: mh | gw"
		  << "\n       Current setting is:  " << chainType << std::endl;
	    return -1;      
	  } else {
	    tcerr << "***Error:  Unrecognized chain type."
		  << "\n       Valid types are: mh | gw"
		  << "\n       Current setting is:  " << chainType << std::endl;
	    return -1;      
	  }
	} else {
	  tcout << "Choose type of chain to run: chain type mh | gw"
		<<"\n      Current setting is:  " << chainType << std::endl;
	}
	break;
      case UNLOAD:
	{
	  if (!nChains) {
	    tcout << "  There are no chains to unload." << std::endl;
	  } else if (nArgs > 2) {
	    if (!prevRange.first)  prevRange.first = 1;
	    if (!prevRange.second) prevRange.second = nChains;
	    RangePair limits = std::make_pair(static_cast<size_t>(1),nChains);
	    StringArray inArgs(nArgs-2);
	    for (size_t i=2; i<nArgs; ++i) inArgs[i-2] = rawArgs[i];
	    IntegerVector removals = 
	      XSparse::getRanges(inArgs, prevRange, limits);
	    chainManager->removeChains(removals);
	    // If the last chain has been removed, 
	    // chainManager will set its isSynchedWithFitParams
	    // flag to false.
	  } else {
	    tcout << "  Range required to specify chain number(s): "
		  << "chain unload <range>" << std::endl;
	  }
	}
	break;
      case WALKERS:
	if (nArgs > 2) {
	  if (!(ssThirdArg >> tmpInput) || !ssThirdArg.eof()) {
	    tcerr << "   Integer walkers parameter required: chain "
		  << "walkers <walkers>" << std::endl;
	    return -1;
	  } else {
	    if (tmpInput % 2 == 1) {
	      string msg("Cannot set chain number of walkers to an odd number.");
	      throw Chain::ChainError(msg);
	    }
	    walkers = tmpInput;
	  }
	} else {
	  tcout << "   Integer walkers parameter required: chain "
		<< "walkers <walkers>" << std::endl;
	}
	break;
      default:
	tcout << "   Subcommand " << subcom << " does not exist" << std::endl;
	XSutility::printValidOptions(tcout, cmd, subComs);
	break; 
      }
    } catch (YellowAlert&) {
      return -1;
    }
  }
  return 0;
}
