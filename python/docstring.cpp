// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Neil Vaytet
#include "docstring.h"


namespace scipp::python {



Docstring::Docstring(const std::string description, const std::string raises, const std::string seealso, const std::string returns, const std::string rtype) : m_description(description), m_raises(raises), m_seealso(seealso), m_returns(returns), m_rtype(rtype) {}


Docstring::Docstring(const std::string description, const std::string raises, const std::string seealso, const std::string returns, const std::string rtype, const std::vector<std::pair<std::string, std::string>> &params) : m_description(description), m_raises(raises), m_seealso(seealso), m_returns(returns), m_rtype(rtype){
  m_params.clear();
  for (auto &p : params)
    m_params.push_back(p);
}

// Docstring::Docstring(const std::string description, const std::string raises, const std::string seealso, const std::string returns, const std::string rtype, const std::pair<std::string, std::string>... params) : m_description(description), m_raises(raises), m_seealso(seealso), m_returns(returns), m_rtype(rtype){
//   m_params.clear();
//   m_params = {params...};
//   // for (auto &p : params)
//   //   m_params.push_back(p);
// }

// void Docstring::set_description(const std::string description);
// void Docstring::set_raises(const std::string raises);
// void Docstring::set_seealso(const std::string seealso);
// void Docstring::set_returns(const std::string returns);
// void Docstring::set_rtype(const std::string rtype);


void Docstring::insert_param(const scipp::index ind, strpair param) {
  m_params.insert(m_params.begin() + ind, param);
}

// void Docstring::erase_param(const scipp::index ind) {
//   m_params.erase(m_params.begin() + ind);
// }

const Docstring Docstring::with_out_arg() {
  Docstring docs(*this);
  docs.m_description += " (in-place)";
  docs.m_rtype += " (View)";
  docs.m_params.push_back({"out", "Output buffer."});
  return docs;
}

const std::string Docstring::to_string() const {

  std::string docstring = m_description + "\n";
  for (auto &p : m_params)
    docstring += ":param " + p.first + ": " + p.second +"\n";
  docstring += ":raises: " + m_raises + "\n";
  docstring += ":seealso: " + m_seealso + "\n";
  docstring += ":return: " + m_returns + "\n";
  docstring += ":rtype: " + m_rtype;
  return docstring;
}


} // namespace scipp::python