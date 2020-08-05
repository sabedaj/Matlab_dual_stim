#include <string>
#include <iostream>

std::string str1 = "Hello=World";

int myInt(std::stoi(str1));
uint16_t myInt16(0);

myInt16 = static_cast<uint16_t>(myInt);

std::cout << myInt16;
