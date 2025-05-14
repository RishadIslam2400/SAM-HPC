#pragma once

#include <stdexcept>
#include <string>
#include <functional>
#include <sstream>

class FailureException : public std::runtime_error
{
public:
    FailureException(const std::string &message) : std::runtime_error(message) {}
};

/* void assertException(const std::string &exceptionClass, const std::function<void()> &callback)
{
    try
    {
        callback();
    }
    catch (const std::exception &e)
    {
        std::string actualClass(typeid(e).name());

        if (actualClass.find(exceptionClass) == std::string::npos)
        {
            throw FailureException("Exception class " + std::string(exceptionClass) + " expected, but " + std::string(actualClass) + " thrown");
        }

        return;
    }

    throw FailureException("Exception expected, but none thrown");
} */

template <typename T>
void assertEquals(const T &a, const T &b, const std::string &message = "")
{
    if (!(a == b))
    {
        std::ostringstream oss;
        if (message.empty())
        {
            oss << "Objects not equal when they should be" << std::endl;
        }
        else
        {
            oss << message << std::endl;
        }

        oss << a << std::endl
            << "expected, but" << std::endl
            << b << " given";

        throw FailureException(oss.str());
    }
}

template <typename X, typename Y>
void assertEquals(const X &a, const Y &b, const std::string &message = "")
{
    if (!(a == b))
    {
        std::ostringstream oss;
        if (message.empty())
        {
            oss << "Objects not equal when they should be" << std::endl;
        }
        else
        {
            oss << message << std::endl;
        }

        oss << a << std::endl
            << "expected, but" << std::endl
            << b << " given";

        throw FailureException(oss.str());
    }
}