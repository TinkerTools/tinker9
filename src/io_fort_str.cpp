#include "io_fort_str.h"
#include "io_text.h"
#include <algorithm>
#include <cassert>
#include <cstring>

TINKER_NAMESPACE_BEGIN
void FortranStringView::copy_with_blank_(char* dst, size_t dstlen,
                                         const char* src, size_t first_n)
{
   if (dst != src) {
      auto m = std::min(dstlen, first_n);
      std::memmove(dst, src, m); // [0, m)
      if (first_n < dstlen) {
         std::fill(&dst[m], &dst[dstlen], ' '); // [m, dstlen)
      }
   }
}

bool FortranStringView::if_eq_(const char* src, size_t len) const
{
   auto lb = this->len_trim();
   auto lc = std::max(lb, len);
   auto buffer = std::string(lc, (char)0);
   // if src is longer, copy b_ to buffer, then compare src and buffer;
   // or copy src to buffer, then compare b_ and buffer
   const char* ptr = b_;
   if (len > lb) {
      copy_with_blank_(&buffer[0], lc, b_, lb);
      ptr = src;
   } else {
      copy_with_blank_(&buffer[0], lc, src, len);
   }
   return !std::strncmp(ptr, buffer.c_str(), lc);
}

FortranStringView::FortranStringView(const char* src, size_t len)
   : b_(const_cast<char*>(src))
   , e_(b_ + len)
{}

FortranStringView::FortranStringView(const char* src)
   : b_(const_cast<char*>(src))
   , e_(b_ + std::strlen(src))
{}

FortranStringView::FortranStringView(const std::string& src)
   : b_(const_cast<char*>(&src[0]))
   , e_(b_ + src.size())
{}

FortranStringView& FortranStringView::operator=(const char* src)
{
   copy_with_blank_(b_, size(), src, std::strlen(src));
   return *this;
}

FortranStringView& FortranStringView::operator=(const std::string& src)
{
   copy_with_blank_(b_, size(), &src[0], src.size());
   return *this;
}

FortranStringView& FortranStringView::operator=(const FortranStringView& src)
{
   copy_with_blank_(b_, size(), src.b_, src.size());
   return *this;
}

bool FortranStringView::operator==(const char* src) const
{
   return if_eq_(src, std::strlen(src));
}

bool FortranStringView::operator==(const std::string& src) const
{
   return if_eq_(src.c_str(), src.size());
}

bool FortranStringView::operator==(const FortranStringView& src) const
{
   return if_eq_(src.b_, src.size());
}

size_t FortranStringView::size() const
{
   return e_ - b_;
}

size_t FortranStringView::len_trim() const
{
   // find the first (char)0
   size_t pos = 0;
   for (; pos < size() && b_[pos] != 0; ++pos)
      ;
   for (; pos > 0 && Text::is_ws(b_[pos - 1]); --pos)
      ;
   return pos;
}

std::string FortranStringView::trim() const
{
   return std::string(b_, b_ + len_trim());
}

FortranStringView FortranStringView::operator()(int begin1, int back1) const
{
   assert(1 <= begin1 && begin1 <= back1 && back1 <= size());
   return FortranStringView(b_ + (begin1 - 1), back1 - begin1 + 1);
}

FortranStringView FortranStringView::operator()(int begin1) const
{
   assert(1 <= begin1 && begin1 <= (e_ - b_));
   return FortranStringView(b_ + (begin1 - 1), e_ - b_);
}
TINKER_NAMESPACE_END
