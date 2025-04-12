#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>
#include <unordered_map>

/// A header-only implementation of parsing command-line options
class CMDOptions {
 private:
  std::string usageMessage;
  std::unordered_map<std::string, std::string> options;

  std::unordered_map<std::string, std::string> allowedOptions;
  std::vector<std::string> allowedOptionsOrder;
  std::unordered_map<std::string, std::string> defaultValues;
  std::unordered_map<std::string, std::vector<std::string> > allowedValues;

  CMDOptions(const CMDOptions&);
  CMDOptions& operator = (const CMDOptions&);
  CMDOptions() {}

 public:
  static std::unique_ptr<CMDOptions> create() {
    return std::unique_ptr<CMDOptions>(new CMDOptions());
  }

  void parse(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
      std::string s(argv[i]);

      if (s == "/?"  || s == "-?" || s == "--help" || s == "-help") {
        usage(argv[0]);
        throw 0;
      }

      parseOption(s);
    }
  }

  void setUsageMessage(const std::string& msg) {
    usageMessage = msg;
  }

  void addAllowedOption(const std::string& optionName, const std::string& defaultValue, const std::string& description) {
    addAllowedOption(optionName, description);
    options[optionName] = defaultValue;
    defaultValues[optionName] = defaultValue;
  }

  void addAllowedOption(const std::string& optionName, const std::string& description) {
    assert(!allowedOptions.count(optionName));
    allowedOptions[optionName] = description;
    allowedOptionsOrder.push_back(optionName);
  }

  void addAllowedValue(const std::string& optionName, const std::string& value) {
    assert(allowedOptions.count(optionName));
    allowedValues[optionName].push_back(value);
  }

  void setStr(const std::string& optionName, const std::string& value) {
    options[optionName] = value;
  }

  std::string getStr(const std::string& optionName) const {
    return getOption(optionName);
  }

  int getInt(const std::string& optionName) const {
    return std::stoi(getOption(optionName));
  }

  void setInt(const std::string& optionName, int value) {
    if (!options.count(optionName)) {
      unrecognizedOption(optionName);
    }

    assert(options.count(optionName));
    setOption(optionName, std::to_string(value));
  }

  bool getBool(const std::string& optionName) const {
    return getOption(optionName) != "false" && getOption(optionName) != "0";
  }

  void setBool(const std::string& optionName, bool value) {
    if (!options.count(optionName)) {
      unrecognizedOption(optionName);
    }

    assert(options.count(optionName));
    setOption(optionName, value ? "true" : "false");
  }

  bool hasOption(const std::string& optionName) const {
    if (!allowedOptions.count(optionName)) {
      unrecognizedOption(optionName);
    }

    return options.count(optionName) > 0;
  }

  bool hasCustomOption(const std::string& optionName) const {
    return options.count(optionName) > 0;
  }

private:
  void parseOption(const std::string& s) {
    const size_t equalIndex = s.find('=');
    std::string name = s.substr(0, equalIndex);

    // custom option
    if (name.substr(0, 2) == "-C") {
      if (name.find('_') != std::string::npos) {
        std::cout << "custom options cannot use underscores: " << name << "\n";
        throw 1;
      }
      name = name.substr(2);
    } else {
      if (!allowedOptions.count(name)) {
        if (equalIndex == std::string::npos && allowedOptions.count("")) {
          options[""] = name;
          return;
        }

        unrecognizedOption(name);
      }
    }

    std::string value = (equalIndex == std::string::npos ? "" : s.substr(equalIndex + 1));
    
    // std::cerr << "parsing option (" << name << ", " << value << ")\n"; 

    if (!options.count(name) || (defaultValues.count(name) && options[name] == defaultValues[name])) {
      options[name] = value;
    }

    if (!allowedValues[name].empty()) {
      bool found = false;
      for (std::string& allowedValue : allowedValues[name]) {
        // check exact match
        if (allowedValue == value) {
          found = true;
          break;
        }
        // check regexp
        if (std::regex_match(value, std::regex(allowedValue))) {
          found = true;
          break;
        }
      }
      if (!found) {
        invalidOption(name);
      }
    }
  }

  std::string getOption(const std::string& optionName) const {
    if (!options.count(optionName)) {
      if (allowedOptions.count(optionName)) {
        unspecifiedOption(optionName);
      }

      unrecognizedOption(optionName);
    }

    assert(options.count(optionName));
    return (*options.find(optionName)).second;
  }

  void setOption(const std::string& optionName, const std::string& value) {
    options[optionName] = value;
  }

  void unspecifiedOption(const std::string& optionName) const {
    std::cout << "required option \"" << optionName << "\" is not specified\n";
    throw 1;
  }

  void unrecognizedOption(const std::string& optionName) const {
    std::cout << "unrecognized option \"" << optionName << "\"\n";
    throw 1;
  }

  void invalidOption(const std::string& optionName) const {
    std::cout << "value \"" << getOption(optionName) << "\" is invalid for option \"" << optionName << "\"\n";
    throw 1;
  }

  void usage(const std::string& program) const {
    if (usageMessage != "") {
      std::cout << usageMessage << "\n";
    } else {
      std::cout << "Usage: " << program << " [options]" << "\n";
    }

    std::cout << "Allowed options:";

    for (auto opt : allowedOptionsOrder) {
      std::string name = allowedOptions.find(opt)->first;

      if (name.length() == 0) {
        continue;
      }

      std::cout << "\n";
      std::cout << "  " << name;

      if (allowedValues.count(name)) {
        auto av = allowedValues.find(name)->second;

        if (!av.empty()) {
          std::cout << "=";
          bool first = true;

          for (std::string s : av)
            if (first)
            { std::cout << "[" << s; first = false; }
            else {
              std::cout << "|" << s;
            }

          std::cout << "]";
        }
      }

      std::cout << "\n";
      std::cout << "  " << allowedOptions.find(opt)->second << "\n";
    }
  }
};
