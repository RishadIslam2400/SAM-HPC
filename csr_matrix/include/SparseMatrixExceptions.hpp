#pragma once

#include <stdexcept>
#include <string>

class InvalidDimensionsException : public std::runtime_error
{
public:
    InvalidDimensionsException(const std::string &message) : std::runtime_error(message) {}
};

class InvalidCoordinatesException : public std::runtime_error
{
public:
    InvalidCoordinatesException(const std::string &message) : std::runtime_error(message) {}
};