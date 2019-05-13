// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Nicolas Senecaut $
// $Authors: Nicolas Senecaut $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_SLIMExport SLIMExport

    @brief Read SLIM FeatureXML and generate a csv fil and a extration of mod 500 peptides for quality control

    @verbinclude UTILS_SLIMExport.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_SLIMExport.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSLIMExport :
  public TOPPBase
{
public:

  TOPPSLIMExport() :
    TOPPBase("SLIMExport", "Read SLIM FeatureXML", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file containing the SLIM featureXML");
    setValidFormats_("in", vector<String>(1, "featureXML"));
    registerOutputFile_("out", "<file>", "", "Output file containing original data from SLIM featureXML");
    setValidFormats_("out", vector<String>(1, "csv"));
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");

    FeatureMap features;
    FeatureXMLFile().load(in, features);
    vector<DataProcessing> dps = features.getDataProcessing();
    bool asymmetric = false;
    for (DataProcessing& dp : dps)
    {
	  if (dp.getMetaValue("parameter: model:type") == "asymmetric")
	  {
		  asymmetric = true;
	  }
	}
    SVOutStream output(out, ",");
	//~ SVOutStream extracted("Extracted_" + out,","); // this file is for quality control of the dynamic feat
	output << "RT"<<"MZ";
	output <<"Charge"<<"Seq"<<"Acession";
	output <<"model height"<<"model_status"<<"model_FWHM"<<"model_center"<<"model_lower"<<"model_upper";
	if (asymmetric)
    {
	    output << "model_EGH_tau";
	    output << "model_EGH_sigma";
    }
    else
    {
	    output << "model_Gauss_sigma";
    }
	output <<"model_area";
	// print Header of isotop
	if (!features[0].getSubordinates().empty()) // check if subordinate is not empty
	{
		int subordinatesSize = features[0].getSubordinates().size();
		for (int i = 0;i<subordinatesSize;i++)
		{
			output << "model height_M" + std::to_string(i)<<"model_status_M" + std::to_string(i);
			output << "model_FWHM_M" + std::to_string(i)<<"model_center_M" + std::to_string(i);
			output << "model_lower_M" + std::to_string(i) << "model_upper_M" + std::to_string(i);
			if (asymmetric)
			{
			    output << "model_EGH_tau_M" + std::to_string(i);
			    output << "model_EGH_sigma_M" + std::to_string(i);
			}
			else
			{
			    output << "model_Gauss_sigma_M" + std::to_string(i);
			}
			output << "model_area_M" + std::to_string(i)<< "Sum_Y_M" + std::to_string(i);
		}	
	}
	output << nl;
	//~ for (unsigned int i = 0;i< features.size();i++)//const Feature& feature : features % 500)
    //~ {
		//~ extracted << "i :" << i << endl;
		//~ extracted << feature.getRT() << feature.getMZ() << nl;
		//~ if (i %500 == 0)
		//~ {
			//~ extracted << "mod i :" << i << nl;
		//~ }
	//~ }
	for (const Feature& feature : features)
    {
		if(feature.getMetaValue("model_status") == "0 (valid)")
		{
			output << feature.getRT() << feature.getMZ();
			//~ for (const auto& pep_id : feature.getPeptideIdentifications()) //Iterate throught PeptideId not executed if empty (because size is 0)
			//~ {
			if (!feature.getPeptideIdentifications().empty())
			{
				const auto& pep_id : feature.getPeptideIdentifications()[0];
				if (!pep_id.getHits().empty())
				{
					const PeptideHit& pep_Hit = pep_id.getHits()[0];
					output << pep_Hit.getCharge(); //Print the charge of the hit
					output << pep_Hit.getSequence(); //Print the sequence of the 1 hit
					if (!pep_Hit.getPeptideEvidences().empty())
					{
						output << pep_Hit.getPeptideEvidences()[0].getProteinAccession(); //Print the sequence of the 1 hit
					}
				}
			}
			output << feature.getMetaValue("model_height");
			output << feature.getMetaValue("model_status");
			output << feature.getMetaValue("model_FWHM");
			output << feature.getMetaValue("model_center");
			output << feature.getMetaValue("model_lower");
			output << feature.getMetaValue("model_upper");
			if (asymmetric)
		    {
			    output << feature.getMetaValue("model_EGH_tau");
			    output << feature.getMetaValue("model_EGH_sigma");
		    }
		    else
		    {
			    output << feature.getMetaValue("model_Gauss_sigma");
		    }
			output << feature.getMetaValue("model_area");
			
			//~ if (!feature.getSubordinates().empty()) // check if subordinate is not empty
			//~ {
				for (const Feature& isotope : feature.getSubordinates()) // loop over the Isotopes
				{
					output << isotope.getMetaValue("model_height");
					output << isotope.getMetaValue("model_status");
					output << isotope.getMetaValue("model_FWHM");
					output << isotope.getMetaValue("model_center");
					output << isotope.getMetaValue("model_lower");
					output << isotope.getMetaValue("model_upper");
					if (asymmetric)
				    {
				    //print correct header
				    output << isotope.getMetaValue("model_EGH_tau");
				    output << isotope.getMetaValue("model_EGH_sigma");
				    }
				    else
				    {
					    //print correct header
					    output << isotope.getMetaValue("model_Gauss_sigma");
				    }
				    output << isotope.getMetaValue("model_area");
				    double sum_y = 0.0;
				    for (const ConvexHull2D& point : isotope.getConvexHulls()) // loop over the YX pairs in the Isotopes
				    {
						//~ sum_y = sum_y + point.getHullPoints()[0];
						for (const auto& line : point.getHullPoints())
						{
							sum_y += line [1];
							//~ output << line[0] << line [1] << nl;
						}
					}
					output << sum_y;
				}
			//~ }
			output << nl;
		}
		
	}

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPSLIMExport tool;
  return tool.main(argc, argv);
}
/// @endcond

