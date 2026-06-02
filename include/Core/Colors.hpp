/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string_view>

namespace gstlrn
{
  struct ColorInfo
  {
    std::string_view name;
    uint32_t hex;
  };

  inline constexpr std::array<ColorInfo, 148> colorTable{
    ColorInfo{"AliceBlue", 0xF0F8FF},
    ColorInfo{"AntiqueWhite", 0xFAEBD7},
    ColorInfo{"Aqua", 0x00FFFF},
    ColorInfo{"Aquamarine", 0x7FFFD4},
    ColorInfo{"Azure", 0xF0FFFF},
    ColorInfo{"Beige", 0xF5F5DC},
    ColorInfo{"Bisque", 0xFFE4C4},
    ColorInfo{"Black", 0x000000},
    ColorInfo{"BlanchedAlmond", 0xFFEBCD},
    ColorInfo{"Blue", 0x0000FF},
    ColorInfo{"BlueViolet", 0x8A2BE2},
    ColorInfo{"Brown", 0xA52A2A},
    ColorInfo{"BurlyWood", 0xDEB887},
    ColorInfo{"CadetBlue", 0x5F9EA0},
    ColorInfo{"Chartreuse", 0x7FFF00},
    ColorInfo{"Chocolate", 0xD2691E},
    ColorInfo{"Coral", 0xFF7F50},
    ColorInfo{"CornflowerBlue", 0x6495ED},
    ColorInfo{"Cornsilk", 0xFFF8DC},
    ColorInfo{"Crimson", 0xDC143C},
    ColorInfo{"Cyan", 0x00FFFF},
    ColorInfo{"DarkBlue", 0x00008B},
    ColorInfo{"DarkCyan", 0x008B8B},
    ColorInfo{"DarkGoldenRod", 0xB8860B},
    ColorInfo{"DarkGray", 0xA9A9A9},
    ColorInfo{"DarkGrey", 0xA9A9A9},
    ColorInfo{"DarkGreen", 0x006400},
    ColorInfo{"DarkKhaki", 0xBDB76B},
    ColorInfo{"DarkMagenta", 0x8B008B},
    ColorInfo{"DarkOliveGreen", 0x556B2F},
    ColorInfo{"DarkOrange", 0xFF8C00},
    ColorInfo{"DarkOrchid", 0x9932CC},
    ColorInfo{"DarkRed", 0x8B0000},
    ColorInfo{"DarkSalmon", 0xE9967A},
    ColorInfo{"DarkSeaGreen", 0x8FBC8F},
    ColorInfo{"DarkSlateBlue", 0x483D8B},
    ColorInfo{"DarkSlateGray", 0x2F4F4F},
    ColorInfo{"DarkSlateGrey", 0x2F4F4F},
    ColorInfo{"DarkTurquoise", 0x00CED1},
    ColorInfo{"DarkViolet", 0x9400D3},
    ColorInfo{"DeepPink", 0xFF1493},
    ColorInfo{"DeepSkyBlue", 0x00BFFF},
    ColorInfo{"DimGray", 0x696969},
    ColorInfo{"DimGrey", 0x696969},
    ColorInfo{"DodgerBlue", 0x1E90FF},
    ColorInfo{"FireBrick", 0xB22222},
    ColorInfo{"FloralWhite", 0xFFFAF0},
    ColorInfo{"ForestGreen", 0x228B22},
    ColorInfo{"Fuchsia", 0xFF00FF},
    ColorInfo{"Gainsboro", 0xDCDCDC},
    ColorInfo{"GhostWhite", 0xF8F8FF},
    ColorInfo{"Gold", 0xFFD700},
    ColorInfo{"GoldenRod", 0xDAA520},
    ColorInfo{"Gray", 0x808080},
    ColorInfo{"Grey", 0x808080},
    ColorInfo{"Green", 0x008000},
    ColorInfo{"GreenYellow", 0xADFF2F},
    ColorInfo{"HoneyDew", 0xF0FFF0},
    ColorInfo{"HotPink", 0xFF69B4},
    ColorInfo{"IndianRed", 0xCD5C5C},
    ColorInfo{"Indigo", 0x4B0082},
    ColorInfo{"Ivory", 0xFFFFF0},
    ColorInfo{"Khaki", 0xF0E68C},
    ColorInfo{"Lavender", 0xE6E6FA},
    ColorInfo{"LavenderBlush", 0xFFF0F5},
    ColorInfo{"LawnGreen", 0x7CFC00},
    ColorInfo{"LemonChiffon", 0xFFFACD},
    ColorInfo{"LightBlue", 0xADD8E6},
    ColorInfo{"LightCoral", 0xF08080},
    ColorInfo{"LightCyan", 0xE0FFFF},
    ColorInfo{"LightGoldenRodYellow", 0xFAFAD2},
    ColorInfo{"LightGray", 0xD3D3D3},
    ColorInfo{"LightGrey", 0xD3D3D3},
    ColorInfo{"LightGreen", 0x90EE90},
    ColorInfo{"LightPink", 0xFFB6C1},
    ColorInfo{"LightSalmon", 0xFFA07A},
    ColorInfo{"LightSeaGreen", 0x20B2AA},
    ColorInfo{"LightSkyBlue", 0x87CEFA},
    ColorInfo{"LightSlateGray", 0x778899},
    ColorInfo{"LightSlateGrey", 0x778899},
    ColorInfo{"LightSteelBlue", 0xB0C4DE},
    ColorInfo{"LightYellow", 0xFFFFE0},
    ColorInfo{"Lime", 0x00FF00},
    ColorInfo{"LimeGreen", 0x32CD32},
    ColorInfo{"Linen", 0xFAF0E6},
    ColorInfo{"Magenta", 0xFF00FF},
    ColorInfo{"Maroon", 0x800000},
    ColorInfo{"MediumAquaMarine", 0x66CDAA},
    ColorInfo{"MediumBlue", 0x0000CD},
    ColorInfo{"MediumOrchid", 0xBA55D3},
    ColorInfo{"MediumPurple", 0x9370DB},
    ColorInfo{"MediumSeaGreen", 0x3CB371},
    ColorInfo{"MediumSlateBlue", 0x7B68EE},
    ColorInfo{"MediumSpringGreen", 0x00FA9A},
    ColorInfo{"MediumTurquoise", 0x48D1CC},
    ColorInfo{"MediumVioletRed", 0xC71585},
    ColorInfo{"MidnightBlue", 0x191970},
    ColorInfo{"MintCream", 0xF5FFFA},
    ColorInfo{"MistyRose", 0xFFE4E1},
    ColorInfo{"Moccasin", 0xFFE4B5},
    ColorInfo{"NavajoWhite", 0xFFDEAD},
    ColorInfo{"Navy", 0x000080},
    ColorInfo{"OldLace", 0xFDF5E6},
    ColorInfo{"Olive", 0x808000},
    ColorInfo{"OliveDrab", 0x6B8E23},
    ColorInfo{"Orange", 0xFFA500},
    ColorInfo{"OrangeRed", 0xFF4500},
    ColorInfo{"Orchid", 0xDA70D6},
    ColorInfo{"PaleGoldenRod", 0xEEE8AA},
    ColorInfo{"PaleGreen", 0x98FB98},
    ColorInfo{"PaleTurquoise", 0xAFEEEE},
    ColorInfo{"PaleVioletRed", 0xDB7093},
    ColorInfo{"PapayaWhip", 0xFFEFD5},
    ColorInfo{"PeachPuff", 0xFFDAB9},
    ColorInfo{"Peru", 0xCD853F},
    ColorInfo{"Pink", 0xFFC0CB},
    ColorInfo{"Plum", 0xDDA0DD},
    ColorInfo{"PowderBlue", 0xB0E0E6},
    ColorInfo{"Purple", 0x800080},
    ColorInfo{"RebeccaPurple", 0x663399},
    ColorInfo{"Red", 0xFF0000},
    ColorInfo{"RosyBrown", 0xBC8F8F},
    ColorInfo{"RoyalBlue", 0x4169E1},
    ColorInfo{"SaddleBrown", 0x8B4513},
    ColorInfo{"Salmon", 0xFA8072},
    ColorInfo{"SandyBrown", 0xF4A460},
    ColorInfo{"SeaGreen", 0x2E8B57},
    ColorInfo{"SeaShell", 0xFFF5EE},
    ColorInfo{"Sienna", 0xA0522D},
    ColorInfo{"Silver", 0xC0C0C0},
    ColorInfo{"SkyBlue", 0x87CEEB},
    ColorInfo{"SlateBlue", 0x6A5ACD},
    ColorInfo{"SlateGray", 0x708090},
    ColorInfo{"SlateGrey", 0x708090},
    ColorInfo{"Snow", 0xFFFAFA},
    ColorInfo{"SpringGreen", 0x00FF7F},
    ColorInfo{"SteelBlue", 0x4682B4},
    ColorInfo{"Tan", 0xD2B48C},
    ColorInfo{"Teal", 0x008080},
    ColorInfo{"Thistle", 0xD8BFD8},
    ColorInfo{"Tomato", 0xFF6347},
    ColorInfo{"Turquoise", 0x40E0D0},
    ColorInfo{"Violet", 0xEE82EE},
    ColorInfo{"Wheat", 0xF5DEB3},
    ColorInfo{"White", 0xFFFFFF},
    ColorInfo{"WhiteSmoke", 0xF5F5F5},
    ColorInfo{"Yellow", 0xFFFF00},
    ColorInfo{"YellowGreen", 0x9ACD32},
  };

  inline bool iequals(std::string_view a, std::string_view b)
  {
    if (a.size() != b.size()) return false;

    for (size_t i = 0; i < a.size(); i++)
    {
      if (std::tolower(static_cast<unsigned char>(a[i]))
          != std::tolower(static_cast<unsigned char>(b[i])))
        return false;
    }
    return true;
  }

  /**
   * @brief Return the Integer code corresponding to the color
   *
   * @param name Color name
   * @return GSTLEARN_EXPORT
   */
  GSTLEARN_EXPORT inline Id nameToId(std::string_view name)
  {
    for (const auto& c: colorTable)
      if (iequals(c.name, name)) return static_cast<Id>(c.hex);

    return static_cast<Id>(-1); // non trouvé
  }

  /**
   * @brief Return the formatted HEXA code corresponding to the color
   *
   * @param name Color name
   * @return GSTLEARN_EXPORT
   */
  GSTLEARN_EXPORT inline String nameToHexa(std::string_view name)
  {
    for (const auto& c: colorTable)
    {
      if (iequals(c.name, name))
      {
        std::ostringstream oss;
        oss << "#" << std::hex << std::uppercase << std::setfill('0')
            << std::setw(6) << c.hex;
        return oss.str();
      }
    }
    return "#INVALID";
  }

  GSTLEARN_EXPORT inline Id nameToPackedRGB(std::string_view name)
  {
    for (const auto& c: colorTable)
      if (iequals(c.name, name)) return static_cast<Id>(c.hex);

    return static_cast<Id>(-1); // couleur non trouvée
  }

  /**
   * @brief Return the Packed RGB code corresponding to an Integer code
   *
   * @param code Integer color code
   * @return Packed RGB code or -1 if not found
   */
  GSTLEARN_EXPORT inline Id intToPackedRGB(Id code)
  {
    for (const auto& c: colorTable)
      if (static_cast<Id>(c.hex) == code) return static_cast<Id>(c.hex);

    return static_cast<Id>(-1);
  }

  /**
   * @brief Return the Packed RGB code corresponding to a HEXA code
   *
   * @param hexa HEXA color code (e.g. "#FFFF00" or "FFFF00")
   * @return Packed RGB code or -1 if invalid
   */
  GSTLEARN_EXPORT inline Id hexaToPackedRGB(const String& hexa)
  {
    std::string value = hexa;

    // Remove optional leading '#'
    if (!value.empty() && value.front() == '#') value.erase(0, 1);

    uint32_t code = 0;

    std::istringstream iss(value);
    iss >> std::hex >> code;

    if (iss.fail()) return static_cast<Id>(-1);

    return static_cast<Id>(code);
  }

  /**
   * @brief Return the Color name corresponding to the Integer code
   *
   * @param code Integer color code
   * @return Color name or an empty string if not found
   */
  GSTLEARN_EXPORT inline String intToName(Id code)
  {
    for (const auto& c: colorTable)
      if (static_cast<Id>(c.hex) == code) return std::string(c.name);

    return "";
  }

  /**
   * @brief Return the Color name corresponding to the HEXA code
   *
   * @param hexa HEXA color code (e.g. "#FFFF00" or "FFFF00")
   * @return Color name or an empty string if not found
   */
  GSTLEARN_EXPORT inline String hexaToName(const String& hexa)
  {
    std::string value = hexa;

    // Remove optional leading '#'
    if (!value.empty() && value.front() == '#') value.erase(0, 1);

    uint32_t code = 0;

    std::istringstream iss(value);
    iss >> std::hex >> code;

    if (iss.fail()) return "";

    for (const auto& c: colorTable)
      if (c.hex == code) return std::string(c.name);

    return "";
  }

  GSTLEARN_EXPORT inline String packedRGBToName(Id code)
  {
    auto rgb = static_cast<uint32_t>(code);

    for (const auto& c: colorTable)
      if (c.hex == rgb) return std::string(c.name);

    return "";
  }
} // namespace gstlrn
