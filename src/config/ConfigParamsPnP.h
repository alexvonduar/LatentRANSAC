#ifndef CONFIGPARAMSPNP_H
#define CONFIGPARAMSPNP_H

#include "ConfigParams.h"

namespace USACConfig
{
    // problem specific/data-related parameters: PnP
    struct PnP
    {
        PnP(): inputFilePath("") // leave blank if not using config file
        {}

        std::string inputFilePath;
    };
}

class ConfigParamsPnP: public ConfigParams
{
public:
    // simple function to read in parameters from config file
    bool initParamsFromConfigFile(std::string& configFilePath);

    USACConfig::PnP pnp;
};

#endif
