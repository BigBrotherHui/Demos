/*
 * File automatically generated by
 * gengen 1.4.1 by Lorenzo Bettini 
 * http://www.gnu.org/software/gengen
 */

#ifndef CHECK_MODES_GEN_CLASS_H
#define CHECK_MODES_GEN_CLASS_H

#include <string>
#include <iostream>

using std::string;
using std::ostream;

class check_modes_gen_class
{
 protected:
  string mode1_given_fields;
  string mode1_name;
  string mode1_options;
  string mode2_given_fields;
  string mode2_name;
  string mode2_options;

 public:
  check_modes_gen_class()
  {
  }
  
  check_modes_gen_class(const string &_mode1_given_fields, const string &_mode1_name, const string &_mode1_options, const string &_mode2_given_fields, const string &_mode2_name, const string &_mode2_options) :
    mode1_given_fields (_mode1_given_fields), mode1_name (_mode1_name), mode1_options (_mode1_options), mode2_given_fields (_mode2_given_fields), mode2_name (_mode2_name), mode2_options (_mode2_options)
  {
  }

  static void
  generate_string(const string &s, ostream &stream, unsigned int indent)
  {
    if (!indent || s.find('\n') == string::npos)
      {
        stream << s;
        return;
      }

    string::size_type pos;
    string::size_type start = 0;
    string ind (indent, ' ');
    while ( (pos=s.find('\n', start)) != string::npos)
      {
        stream << s.substr (start, (pos+1)-start);
        start = pos+1;
        if (start+1 <= s.size ())
          stream << ind;
      }
    if (start+1 <= s.size ())
      stream << s.substr (start);
  }

  void set_mode1_given_fields(const string &_mode1_given_fields)
  {
    mode1_given_fields = _mode1_given_fields;
  }

  void set_mode1_name(const string &_mode1_name)
  {
    mode1_name = _mode1_name;
  }

  void set_mode1_options(const string &_mode1_options)
  {
    mode1_options = _mode1_options;
  }

  void set_mode2_given_fields(const string &_mode2_given_fields)
  {
    mode2_given_fields = _mode2_given_fields;
  }

  void set_mode2_name(const string &_mode2_name)
  {
    mode2_name = _mode2_name;
  }

  void set_mode2_options(const string &_mode2_options)
  {
    mode2_options = _mode2_options;
  }

  void generate_check_modes(ostream &stream, unsigned int indent = 0);
  
};

#endif // CHECK_MODES_GEN_CLASS_H
