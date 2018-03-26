break_into_subregions = function(x,y, m_y, m_x,...){
  y_range = diff(range(y))
  x_range = diff(range(x))
  x_subrange = x_range/m_x
  y_subrange = y_range/m_y
  x_breaks = min(x)+  c(0:m_x)*x_subrange
  y_breaks = min(y) + c(0:m_y)*y_subrange
  x_bin = cut(x, x_breaks, labels=c(1:m_x), include.lowest=T)
  y_bin = cut(y, y_breaks, labels=c(1:m_y), include.lowest=T)
  bin = paste0(x_bin, y_bin)
  ;
  data.frame(x,y,bin)
}