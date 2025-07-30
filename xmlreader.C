#include "TXMLEngine.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

void FindNodeInfo(TXMLEngine &xml, XMLNodePointer_t node, std::vector<std::string> &enabledChannelsVector, std::map<std::string, std::map<std::string, std::string>> &channelInfoMap) {
   const char *nodeName = xml.GetNodeName(node);
   const char *nodeContent = xml.GetNodeContent(node);

   if (std::string(nodeName) == "channel") {
      XMLNodePointer_t index = xml.GetChild(node);
      if (index != nullptr) {
         std::map<std::string, std::string> currentChannelInfo;
         std::string indexNumber = std::string(xml.GetNodeContent(index));
         // allChannels.push_back(indexNumber);
         XMLNodePointer_t entry = xml.GetChild(xml.GetNext(index));
         while (entry != nullptr) { // Get all entry nodes of the channel
            XMLNodePointer_t key = xml.GetChild(entry);
            XMLNodePointer_t value = xml.GetNext(key);
            if (key != nullptr && value != nullptr) { // There are key nodes with no value nodes
               currentChannelInfo[std::string(xml.GetNodeContent(key))] = std::string(xml.GetNodeContent(value));
               if (std::string(xml.GetNodeContent(key)) == "SRV_PARAM_CH_ENABLED" && std::string(xml.GetNodeContent(value)) == "true") // Checks if the channel is enabled
                  enabledChannelsVector.push_back(indexNumber);
               channelInfoMap[indexNumber] = currentChannelInfo;
            }
            entry = xml.GetNext(entry);
         }
      }
   }

   XMLNodePointer_t child = xml.GetChild(node);
   while (child != nullptr) {
      FindNodeInfo(xml, child, enabledChannelsVector, channelInfoMap);
      child = xml.GetNext(child);
   }
}

void xmlreader(const char *filePath = "file/path/to/settings.xml") {
   // Create engine
   TXMLEngine xml;

   // Parse the settings.xml file
   XMLDocPointer_t xmldoc = xml.ParseFile(filePath);
   if (!xmldoc)
      std::cerr << "File not parsed" << std::endl;

   // Find root node
   XMLNodePointer_t rootNode = xml.DocGetRootElement(xmldoc);

   // std::vector<std::string> allChannels;
   std::vector<std::string> enabledChannels;
   enabledChannels.clear(); // in case the vector isn't empty
   std::map<std::string, std::map<std::string, std::string>> channelInfo;
   channelInfo.clear();

   // Get channel info
   FindNodeInfo(xml, rootNode, enabledChannels, channelInfo);

   // Filtering out channels that are not enabled
   std::map<std::string, std::map<std::string, std::string>> filteredChannelInfo;

   for (const auto &outerPair : channelInfo) {
      const std::string &channelIndex = outerPair.first;

      if (std::find(enabledChannels.begin(), enabledChannels.end(), channelIndex) != enabledChannels.end()) {
         filteredChannelInfo[channelIndex] = outerPair.second;
      }
   }
   channelInfo = filteredChannelInfo;

   // Prints out the names of all enabled channels
   std::cout << "The " << enabledChannels.size() << " enabled channels are: ";
   for (auto &index : enabledChannels) {
      std::cout << "channel " << index << " (" << channelInfo[index]["SW_PARAMETER_CH_LABEL"] << ") ";
   }
   std::cout << std::endl;

   // Setting calibration parameters
   for (auto &index : enabledChannels) {
      if (channelInfo[index]["SW_PARAMETER_CH_LABEL"] == "Up Half") {
         float gainUp = std::stof(channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_P0"]);
         float offsetUp = std::stof(channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_P1"]);
         float quadraticUp = std::stof(channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_P2"]);
         std::string units = channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_UDM"];
         std::cout << "For Up Half:" << std::endl;
         std::cout << "Gain: " << gainUp << ", Offset: " << offsetUp << ", Quadratic: " << quadraticUp << ", Units: " << units << std::endl;
      } else if (channelInfo[index]["SW_PARAMETER_CH_LABEL"] == "Down Half") {
         float gainDown = std::stof(channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_P0"]);
         float offsetDown = std::stof(channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_P1"]);
         float quadraticDown = std::stof(channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_P2"]);
         std::string units = channelInfo[index]["SW_PARAMETER_CH_ENERGY_CALIBRATION_UDM"];
         std::cout << "For Down Half:" << std::endl;
         std::cout << "Gain: " << gainDown << ", Offset: " << offsetDown << ", Quadratic: " << quadraticDown << ", Units: " << units << std::endl;
      }
   }

   // Release memory before exit
   xml.FreeDoc(xmldoc);
   // Delete engine
   xml.Delete();
}
