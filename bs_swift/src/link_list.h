#ifndef __LINK_LIST_H_
#define __LINK_LIST_H_
/*!
 * \file link_list.h
 * \brief double linked list declaration
 * \author Khait Mark
 * \date 2006-07-24
 */

#define LIST_HEAD -1
#define LIST_TAIL -2

namespace blue_sky
{
/**
 * @brief double linked list for AMG
 */
struct double_linked_list
  {
    t_long                    data;         //!< data
    t_long                    head;
    t_long                    tail;
    struct double_linked_list *next_elt;     //!< pointer to the next node
    struct double_linked_list *prev_elt;     //!< pointer to the previous node
  };



/**
 * @brief remove node
 *
 * @param element_ptr -- node
 */
void dispose_elt (double_linked_list *element_ptr )
{
  delete element_ptr;
}



/**
 * @brief remove_point:   removes a point from the lists
 *
 * @param LoL_head_ptr
 * @param LoL_tail_ptr
 * @param meassure
 * @param index
 * @param lists
 * @param where
 */
void remove_point (double_linked_list **LoL_head_ptr, 
                   double_linked_list **LoL_tail_ptr,
                   t_long meassure, 
                   t_long index, 
                   t_long *lists, 
                   t_long *where)
{
  double_linked_list *LoL_head = *LoL_head_ptr;
  double_linked_list *LoL_tail = *LoL_tail_ptr;
  double_linked_list *list_ptr;

  list_ptr =  LoL_head;


  do
    {
      if (meassure == list_ptr->data)
        {

          /* point to be removed is only point on list,
             which must be destroyed */
          if (list_ptr->head == index && list_ptr->tail == index)
            {
              /* removing only list, so num_left better be 0! */
              if (list_ptr == LoL_head && list_ptr == LoL_tail)
                {
                  LoL_head = 0;
                  LoL_tail = 0;
                  dispose_elt(list_ptr);

                  *LoL_head_ptr = LoL_head;
                  *LoL_tail_ptr = LoL_tail;
                  return;
                }
              else if (LoL_head == list_ptr) /*removing 1st (max_meassure) list */
                {
                  list_ptr -> next_elt -> prev_elt = 0;
                  LoL_head = list_ptr->next_elt;
                  dispose_elt(list_ptr);

                  *LoL_head_ptr = LoL_head;
                  *LoL_tail_ptr = LoL_tail;
                  return;
                }
              else if (LoL_tail == list_ptr)     /* removing last list */
                {
                  list_ptr -> prev_elt -> next_elt = 0;
                  LoL_tail = list_ptr->prev_elt;
                  dispose_elt(list_ptr);

                  *LoL_head_ptr = LoL_head;
                  *LoL_tail_ptr = LoL_tail;
                  return;
                }
              else
                {
                  list_ptr -> next_elt -> prev_elt = list_ptr -> prev_elt;
                  list_ptr -> prev_elt -> next_elt = list_ptr -> next_elt;
                  dispose_elt(list_ptr);

                  *LoL_head_ptr = LoL_head;
                  *LoL_tail_ptr = LoL_tail;
                  return;
                }
            }
          else if (list_ptr->head == index)      /* index is head of list */
            {
              list_ptr->head = lists[index];
              where[lists[index]] = LIST_HEAD;
              return;
            }
          else if (list_ptr->tail == index)      /* index is tail of list */
            {
              list_ptr->tail = where[index];
              lists[where[index]] = LIST_TAIL;
              return;
            }
          else                              /* index is in middle of list */
            {
              lists[where[index]] = lists[index];
              where[lists[index]] = where[index];
              return;
            }
        }
      list_ptr = list_ptr -> next_elt;
    }
  while (list_ptr != 0);

  return;
}

/**
 * @brief create_elt() : Create an element using Item for its data field
 *
 * @param Item
 *
 * @return node
 */
double_linked_list *create_elt (t_long Item)
{
  double_linked_list *new_elt_ptr;


  /* Allocate memory space for the new node.
   * return with error if no space available
   */

  if ( (new_elt_ptr = (new double_linked_list)) == 0)
    {
      return 0;
    }
  else

    /*   new_elt_ptr = hypre_CTAlloc(hypre_LinkList, 1); */

    {
      new_elt_ptr -> data = Item;
      new_elt_ptr -> next_elt = 0;
      new_elt_ptr -> prev_elt = 0;
      new_elt_ptr -> head = LIST_TAIL;
      new_elt_ptr -> tail = LIST_HEAD;
    }

  return (new_elt_ptr);
}

/**
 * @brief enter_on_lists   places point in new list
 *
 * @param LoL_head_ptr
 * @param LoL_tail_ptr
 * @param meassure
 * @param index
 * @param lists
 * @param where
 */
void enter_on_lists (double_linked_list **LoL_head_ptr, 
                   double_linked_list **LoL_tail_ptr, 
                   t_long meassure, 
                   t_long index, 
                   t_long *lists, 
                   t_long *where)
{
  double_linked_list *LoL_head = *LoL_head_ptr;
  double_linked_list *LoL_tail = *LoL_tail_ptr;

  double_linked_list *list_ptr;
  double_linked_list *new_ptr;

  t_long         old_tail;

  list_ptr =  LoL_head;

  if (LoL_head == 0)   /* no lists exist yet */
    {
      new_ptr = create_elt (meassure);
      new_ptr->head = index;
      new_ptr->tail = index;
      lists[index] = LIST_TAIL;
      where[index] = LIST_HEAD;
      LoL_head = new_ptr;
      LoL_tail = new_ptr;

      *LoL_head_ptr = LoL_head;
      *LoL_tail_ptr = LoL_tail;
      return;
    }
  else
    {
      do
        {
          if (meassure > list_ptr->data)
            {
              new_ptr = create_elt (meassure);
              new_ptr->head = index;
              new_ptr->tail = index;
              lists[index] = LIST_TAIL;
              where[index] = LIST_HEAD;

              if ( list_ptr->prev_elt != 0)
                {
                  new_ptr->prev_elt            = list_ptr->prev_elt;
                  list_ptr->prev_elt->next_elt = new_ptr;
                  list_ptr->prev_elt           = new_ptr;
                  new_ptr->next_elt            = list_ptr;
                }
              else
                {
                  new_ptr->next_elt  = list_ptr;
                  list_ptr->prev_elt = new_ptr;
                  new_ptr->prev_elt  = 0;
                  LoL_head = new_ptr;
                }

              *LoL_head_ptr = LoL_head;
              *LoL_tail_ptr = LoL_tail;
              return;
            }
          else if (meassure == list_ptr->data)
            {
              old_tail = list_ptr->tail;
              lists[old_tail] = index;
              where[index] = old_tail;
              lists[index] = LIST_TAIL;
              list_ptr->tail = index;
              return;
            }

          list_ptr = list_ptr->next_elt;
        }
      while (list_ptr != 0);

      new_ptr = create_elt (meassure);
      new_ptr->head = index;
      new_ptr->tail = index;
      lists[index] = LIST_TAIL;
      where[index] = LIST_HEAD;
      LoL_tail->next_elt = new_ptr;
      new_ptr->prev_elt = LoL_tail;
      new_ptr->next_elt = 0;
      LoL_tail = new_ptr;

      *LoL_head_ptr = LoL_head;
      *LoL_tail_ptr = LoL_tail;
      return;
    }
}

}
#endif //__LINK_LIST_H_
