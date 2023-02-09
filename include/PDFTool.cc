#include "PDFTool.hh"

std::pair <int, int> PDF_ID_range;

std::pair< LHAPDF::PDF*, LHAPDF::PDF* > alphaS_sets;
std::pair <int, int> alphaS_IDs;
LHAPDF::PDF const *pdfProd;

rvec<float> calc_pdf_weights(   const float& q2scale,
                                const float& x1,
                                const float& x2,
                                const int& id1,
                                const int& id2) {
    rvec<float> pdf_weights;
    
    float prodWeight = pdfProd->xfxQ( id1, x1, q2scale ) * pdfProd->xfxQ( id2, x2, q2scale );
    
    for( auto& PDFSet : config::init_pdf_sets ) {
        pdf_weights.push_back( PDFSet->xfxQ( id1, x1, q2scale ) * PDFSet->xfxQ( id2, x2, q2scale ) / prodWeight );
    }
    
    return pdf_weights;
}

std::pair<float, float> calc_as_weights(const float& q2scale,
                                        const float& x1,
                                        const float& x2,
                                        const int& id1,
                                        const int& id2) {
                                                
    
    float prodWeight = pdfProd->xfxQ( id1, x1, q2scale ) * pdfProd->xfxQ( id2, x2, q2scale );
    
    float as_weight_down = alphaS_sets.first->xfxQ( id1, x1, q2scale ) * alphaS_sets.first->xfxQ( id2, x2, q2scale ) / prodWeight;
    float as_weight_up = alphaS_sets.second->xfxQ( id1, x1, q2scale ) * alphaS_sets.second->xfxQ( id2, x2, q2scale ) / prodWeight;

    return std::make_pair( as_weight_down, as_weight_up);
}

RNode init_PDFs(RNode& df) {
    // check if info is already present
    std::vector<std::string> colNames = df.GetColumnNames();
    if (std::find(colNames.begin(), colNames.end(), "LHEPdfWeight_def") == colNames.end()) {
		
        ///////////////////////////////////////////////////////////////////////////
        // this is initialization of pdf, which only needs to be done once
        if (!config::pdf_is_initialized) {
            if (config::pdf_set_name != "") {
                PDF_ID_range = std::make_pair( config::pdf_setid, (config::pdf_setid + config::pdf_nweights - 1) );
                alphaS_IDs = std::make_pair( 0, 0 );
            }
            else {
                PDF_ID_range = std::make_pair( 0, 0 );
                alphaS_IDs = std::make_pair( 0, 0 );
            }

            for ( int setID = PDF_ID_range.first; setID <= PDF_ID_range.second; setID++ ) {
                config::init_pdf_sets.push_back( LHAPDF::mkPDF( setID ) );
            }
            
            if( alphaS_IDs.first and alphaS_IDs.second) {
                alphaS_sets = std::make_pair( LHAPDF::mkPDF( alphaS_IDs.first ),
                                              LHAPDF::mkPDF( alphaS_IDs.second ) );
            }
            
            pdfProd = LHAPDF::mkPDF( config::pdf_prod_set_name, 0 );
            config::pdf_is_initialized = true;
        }
        //////////////////////////////////////////////////////////////////////////
        
        
        // write info and length to LHEPdfWeight columns
        
        df = df.Define("LHEPdfWeight_def", calc_pdf_weights, {"Generator_scalePDF", "Generator_x1", "Generator_x2", "Generator_id1", "Generator_id2"});
        if( alphaS_IDs.first and alphaS_IDs.second) {
            df = df.Define("alphaS_weights", calc_as_weights, {"Generator_scalePDF", "Generator_x1", "Generator_x2", "Generator_id1", "Generator_id2"});
        }
        df = df.Define("nLHEPdfWeight_def", [](rvec<float> pdf_weights){return pdf_weights.size();}, {"LHEPdfWeight_def"});
	}
    return df;
}

