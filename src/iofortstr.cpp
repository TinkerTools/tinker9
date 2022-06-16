#include "tool/iofortstr.h"
#include "tool/iotext.h"
#include <algorithm>

namespace tinker {
void FortranStringView::copyWithBlank(char* dst, size_t dstlen, const char* src, size_t first_n)
{
   if (dst != src) {
      auto m = std::min(dstlen, first_n);
      std::memmove(dst, src, m); // [0, m)
      if (first_n < dstlen) {
         std::fill(&dst[m], &dst[dstlen], ' '); // [m, dstlen)
      }
   }
}

bool FortranStringView::ifEq(const char* src, size_t len) const
{
   auto lb = this->len_trim();
   auto lc = std::max(lb, len);
   auto buffer = std::string(lc, (char)0);
   // If src is longer, copy m_b to buffer, then compare src and buffer;
   // or copy src to buffer, then compare m_b and buffer.
   const char* ptr = m_b;
   if (len > lb) {
      copyWithBlank(&buffer[0], lc, m_b, lb);
      ptr = src;
   } else {
      copyWithBlank(&buffer[0], lc, src, len);
   }
   return 0 == std::strncmp(ptr, buffer.c_str(), lc);
}

size_t FortranStringView::size() const
{
   return m_e - m_b;
}

FortranStringView::FortranStringView(char* src, size_t len)
   : m_b(src)
   , m_e(m_b + len)
{}

FortranStringView& FortranStringView::operator=(const char* src)
{
   copyWithBlank(m_b, size(), src, std::strlen(src));
   return *this;
}

FortranStringView& FortranStringView::operator=(const std::string& src)
{
   copyWithBlank(m_b, size(), &src[0], src.size());
   return *this;
}

bool FortranStringView::operator==(const char* src) const
{
   return ifEq(src, std::strlen(src));
}

bool FortranStringView::operator==(const std::string& src) const
{
   return ifEq(src.c_str(), src.size());
}

bool FortranStringView::operator==(const FortranStringView& src) const
{
   return ifEq(src.m_b, src.size());
}

size_t FortranStringView::len_trim() const
{
   // Find the first (char)0.
   size_t pos = 0;
   for (; pos < size() && m_b[pos] != 0; ++pos)
      ;
   for (; pos > 0 && Text::isWhiteSpace(m_b[pos - 1]); --pos)
      ;
   return pos;
}

std::string FortranStringView::trim() const
{
   return std::string(m_b, m_b + len_trim());
}

FortranStringView FortranStringView::operator()(int begin1, int back1) const
{
   return FortranStringView(m_b + (begin1 - 1), back1 - begin1 + 1);
}

FortranStringView FortranStringView::operator()(int begin1) const
{
   return FortranStringView(m_b + (begin1 - 1), m_e - m_b);
}
}
