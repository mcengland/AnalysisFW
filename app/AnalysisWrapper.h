
/**
 * @brief Wrapper classes are defined here to load them in python
 * Wrapper cannot return references or custom classes.
 */
#include <memory>
#include <TTree.h>
#include <string>
#include "CLoop.h"


class CLoopWrapper {
    public:
    explicit CLoopWrapper(TTree* tree, const std::string& sampleName) {
        m_cloop = std::make_shared<CLoop>(tree, sampleName);
    }

    void Loop(double lumFactor, int z_sample, std::string key){
        m_cloop->Loop(lumFactor, z_sample, key);
    }
    
    private:
    std::shared_ptr<CLoop> m_cloop;

};
