#ifndef GEMDigi_GEMDigiPreReco_h
#define GEMDigi_GEMDigiPreReco_h

/** \class GEMDigiPreReco
 *
 * Digi for GEM
 *  
 * \author Marcello Maggi
 *
 */

#include <boost/cstdint.hpp>
#include <iosfwd>

class GEMDigiPreReco{

public:
//  explicit GEMDigiPreReco (float x, float y, float ex, float ey, float corr, float tof);
  explicit GEMDigiPreReco (float x, float y, float ex, float ey, float corr, float tof, int pdgid);
  GEMDigiPreReco ();

  bool operator==(const GEMDigiPreReco& digi) const;
  bool operator!=(const GEMDigiPreReco& digi) const;
  bool operator<(const GEMDigiPreReco& digi) const;

  float x() const { return x_; }
  float y() const { return y_; }
  float ex() const { return ex_; }
  float ey() const { return ey_; }
  float corr() const { return corr_; }
  float tof() const { return tof_;}
//cesare changes
  int pdgid() const { return pdgid_;}
  void print() const;

private:
  float x_;
  float y_;
  float ex_;
  float ey_;
  float corr_;
  float tof_;
//cesare changes
  int pdgid_;
};

std::ostream & operator<<(std::ostream & o, const GEMDigiPreReco& digi);

#endif

