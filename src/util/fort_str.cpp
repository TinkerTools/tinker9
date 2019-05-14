#include "util/fort_str.h"
#include <algorithm>

TINKER_NAMESPACE_BEGIN
namespace detail_ {
static const char BLANK = ' ';
}

void FortranStringView::copy_with_blank_(char* _dst, size_t _dstlen,
                                         const char* _src, size_t _first_n) {
  if (_dst != _src) {
    auto m = std::min(_dstlen, _first_n);
    std::memmove(_dst, _src, m); // [0, m)
    if (_first_n < _dstlen) {
      std::fill(&_dst[m], &_dst[_dstlen], detail_::BLANK); // [m, _dstlen)
    }
  }
}

bool FortranStringView::if_eq_(const char* _src, size_t _len) const {
  auto lb = this->len_trim();
  auto lc = std::max(lb, _len);
  auto buffer = std::string(lc, (char)0);
  // If _src is longer, copy b_ to buffer, then compare _src and buffer;
  // or copy _src to buffer, then compare b_ and buffer.
  const char* ptr = b_;
  if (_len > lb) {
    copy_with_blank_(&buffer[0], lc, b_, lb);
    ptr = _src;
  } else {
    copy_with_blank_(&buffer[0], lc, _src, _len);
  }
  return !std::strncmp(ptr, buffer.c_str(), lc);
}

//======================================================================

FortranStringView::FortranStringView(const char* _src, size_t _len)
    : b_(const_cast<char*>(_src)), e_(b_ + _len) {}

FortranStringView::FortranStringView(const char* _src)
    : b_(const_cast<char*>(_src)), e_(b_ + std::strlen(_src)) {}

FortranStringView::FortranStringView(const std::string& _src)
    : b_(const_cast<char*>(&_src[0])), e_(b_ + _src.size()) {}

//======================================================================

FortranStringView& FortranStringView::operator=(const char* _src) {
  copy_with_blank_(b_, size(), _src, std::strlen(_src));
  return *this;
}

FortranStringView& FortranStringView::operator=(const std::string& _src) {
  copy_with_blank_(b_, size(), &_src[0], _src.size());
  return *this;
}

FortranStringView& FortranStringView::operator=(const FortranStringView& _src) {
  copy_with_blank_(b_, size(), _src.b_, _src.size());
  return *this;
}

//======================================================================

bool FortranStringView::operator==(const char* __src) const {
  return if_eq_(__src, std::strlen(__src));
}

bool FortranStringView::operator==(const std::string& __src) const {
  return if_eq_(__src.c_str(), __src.size());
}

bool FortranStringView::operator==(const FortranStringView& __src) const {
  return if_eq_(__src.b_, __src.size());
}

//======================================================================

size_t FortranStringView::size() const { return e_ - b_; }

size_t FortranStringView::len_trim() const {
  // find the first (char)0
  size_t pos = 0;
  for (; pos < size() && b_[pos] != 0; ++pos)
    ;
  for (; pos > 0 && b_[pos - 1] == detail_::BLANK; --pos)
    ;
  return pos;
}

std::string FortranStringView::trim() const {
  return std::string(b_, b_ + len_trim());
}

FortranStringView FortranStringView::operator()(int __begin1,
                                                int __back1) const {
  assert(1 <= __begin1 && __begin1 <= __back1 && __back1 <= size());
  return FortranStringView(b_ + (__begin1 - 1), __back1 - __begin1 + 1);
}

TINKER_NAMESPACE_END
